// rc_c152_sim.cpp
// 1:6 scale Cessna 152 RC model: takeoff / climb / cruise simulator + design optimizer.
//
// Build:  g++ -O2 -std=c++17 rc_c152_sim.cpp -o rc_sim
// Run:    ./rc_sim
//
// MODEL (3-DOF longitudinal, SI units):
//   states: x, z (height), vx, vz, pitch angle theta, pitch rate q
//   forces: lift, drag (parabolic polar + stall + ground effect), prop thrust
//           (falls off linearly with airspeed), weight, gear rolling friction
//   wind:   steady wind (power-law shear with height) at any angle to the runway
//           (headwind component acts in the sim, crosswind component is reported as
//           required crab angle) + Gauss-Markov turbulence in x and z
//   pilot:  simple autopilot: full-throttle roll, rotate at 1.15*Vstall, climb at
//           1.4*Vstall, then altitude/airspeed hold. Pitch loop is gain-scheduled.
//
// LIMITS: no lateral-directional dynamics (no roll/yaw), no flaps, tail lift is only
// used in the pitching moment. Treat results as a design-stage comparison tool, not
// a flight-test replacement.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double G = 9.80665;
constexpr double PI = 3.14159265358979323846;
inline double d2r(double d) { return d * PI / 180.0; }
inline double r2d(double r) { return r * 180.0 / PI; }
inline double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// ------------------------------------------------------------------ data types
struct Design {
    double mass = 2.5;          // kg  TOTAL (worst case). 1 lb payload = 0.454 kg of it
    double payload = 0.454;     // kg  (info only)
    double span = 1.697;        // m   (C-152: 10.17 m / 6)
    double chord = 0.2447;      // m   (rectangular wing, C-152: 1.47 m / 6)
    double length = 1.223;      // m   fuselage length (7.34 m / 6)
    double tailArea = 0.0556;   // m^2 horizontal tail (approx 2.0 m^2 / 36)
    double tailArm = 0.73;      // m   CG to tail aerodynamic centre
    double cgFrac = 0.28;       // CG position, fraction of chord behind wing leading edge
    double thrustStatic = 9.0;  // N   static thrust of motor/prop (~0.6 x weight)
    double vPitch = 22.0;       // m/s prop "pitch speed" (thrust -> 0 at this airspeed)
    double cd0 = 0.045;         // zero-lift drag (fixed gear, struts, low Re)
    double clMax = 1.30;        // max lift coefficient (low Re, no flaps)
    double alphaZeroDeg = -3.5; // zero-lift AoA vs fuselage line (airfoil -2 deg, wing incidence +1.5)
    double cmac = -0.045;       // wing pitching moment coeff about AC (cambered airfoil)
    double oswald = 0.80;       // span efficiency
};

struct Env {
    double windSpeed = 4.0;     // m/s steady wind at 10 m height
    double windAngleDeg = 0.0;  // 0 = direct headwind, 90 = from the right (pure crosswind)
    double gustFrac = 0.30;     // turbulence sigma as fraction of wind speed
    double rho = 1.225;         // kg/m^3
    double runwayLen = 30.0;    // m available for takeoff roll
    double rollFriction = 0.04; // paved ~0.03-0.04, short grass ~0.06-0.10
    double cruiseAlt = 25.0;    // m target altitude
    double nLimit = 5.0;        // structural load-factor limit
};

struct Aero {
    double S, AR, CLa, CLat, Vh, Cma, Cm0, Cmde, Cmq, k, alpha0, alphaStall, Iyy, Vstall, SM, Vcruise;
};

Aero computeAero(const Design& d, double rho) {
    Aero a;
    a.S = d.span * d.chord;
    a.AR = d.span * d.span / a.S;
    const double a0 = 2.0 * PI;
    a.CLa = a0 / (1.0 + a0 / (PI * d.oswald * a.AR));
    a.CLat = 3.5;                                    // low-AR tail lift slope per rad
    a.Vh = d.tailArea * d.tailArm / (a.S * d.chord); // tail volume coefficient
    const double deda = 2.0 * a.CLa / (PI * a.AR);   // downwash gradient
    const double eta = 0.9, tau = 0.5, h = d.cgFrac;
    a.alpha0 = d2r(d.alphaZeroDeg);
    a.Cma = a.CLa * (h - 0.25) - eta * a.Vh * a.CLat * (1.0 - deda);
    a.Cm0 = d.cmac - a.CLa * a.alpha0 * (h - 0.25) - eta * a.Vh * a.CLat * deda * a.alpha0;
    a.Cmde = -eta * a.Vh * a.CLat * tau;
    a.Cmq = -2.0 * eta * a.CLat * a.Vh * (d.tailArm / d.chord);
    a.k = 1.0 / (PI * d.oswald * a.AR);
    a.alphaStall = a.alpha0 + d.clMax / a.CLa;
    const double ky = 0.25 * d.length;
    a.Iyy = d.mass * ky * ky;
    a.Vstall = std::sqrt(2.0 * d.mass * G / (rho * a.S * d.clMax));
    a.SM = -a.Cma / a.CLa; // static margin (fraction of chord); want ~0.05-0.20
    a.Vcruise = 1.8 * a.Vstall;
    return a;
}

double liftCoeff(double alpha, const Design& d, const Aero& a, bool& stalled) {
    const double da = alpha - a.alpha0, dS = d.clMax / a.CLa;
    stalled = std::fabs(da) > dS;
    if (!stalled) return a.CLa * da;
    const double s = da > 0 ? 1.0 : -1.0;
    return s * std::max(0.55 * d.clMax, d.clMax - 1.8 * (std::fabs(da) - dS));
}

struct SimResult {
    bool ok = false;
    std::string reason = "";
    bool lifted = false, reachedCruise = false;
    double rollDist = 0, liftoffV = 0, liftoffT = 0, tToAlt = 0;
    double altRMS = 99, maxN = 0, minAlphaMargin = 9, maxAlt = 0, crabDeg = 0, LDcruise = 0;
    struct Row { double t, x, z, V, alphaDeg, thetaDeg, n, thr, de; };
    std::vector<Row> trace;
};

enum Phase { ROLL, ROTATE, CLIMB, CRUISE };

// ------------------------------------------------------------------ simulation
SimResult simulate(const Design& d, const Env& e, unsigned seed, bool wantTrace = false) {
    SimResult R;
    const Aero a = computeAero(d, e.rho);
    const double dt = 0.005, Tmax = 50.0;
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);

    const double hw = e.windSpeed * std::cos(d2r(e.windAngleDeg));
    const double cw = e.windSpeed * std::sin(d2r(e.windAngleDeg));
    const double sigma = e.gustFrac * e.windSpeed;
    const double tauG = 0.8;

    double x = 0, z = 0, vx = 0, vz = 0, th = 0, q = 0;
    double gx = 0, gz = 0, intTh = 0, intAlt = 0;
    Phase ph = ROLL;
    double tLift = -1, stallRun = 0, sumErr2 = 0, tCruise = -1;
    long nErr = 0;
    double nextTrace = 0.0;
    const double Vr = 1.15 * a.Vstall, Vclimb = 1.4 * a.Vstall;
    const double rotTheta = d2r(9.0), deMax = d2r(25.0);

    {   // cruise L/D (info)
        const double CL = 2 * d.mass * G / (e.rho * a.S * a.Vcruise * a.Vcruise);
        R.LDcruise = CL / (d.cd0 + a.k * CL * CL);
    }

    for (long step = 0; step * dt < Tmax; ++step) {
        const double t = step * dt;

        // ---- wind
        gx += -gx / tauG * dt + sigma * std::sqrt(2 * dt / tauG) * nd(rng);
        gz += -gz / tauG * dt + 0.5 * sigma * std::sqrt(2 * dt / tauG) * nd(rng);
        const double prof = std::pow(std::max(z, 0.3) / 10.0, 0.143);
        const double wx = -hw * prof + gx;           // wind blowing toward -x = headwind
        const double wz = (z > 0.05) ? gz : 0.0;
        const double Vrx = vx - wx, Vrz = vz - wz;
        const double V = std::max(0.1, std::hypot(Vrx, Vrz));
        const double gam = std::atan2(Vrz, Vrx);
        const double alpha = th - gam;

        // ---- aerodynamics
        bool stalled = false;
        const double CL = liftCoeff(alpha, d, a, stalled);
        double CD = d.cd0 + a.k * CL * CL;
        const double hAbove = z + 0.15;
        const double ge = std::pow(16 * hAbove / d.span, 2);
        CD = d.cd0 + a.k * CL * CL * (ge / (1 + ge)); // ground effect on induced drag
        if (stalled) CD += 2.0 * (std::fabs(alpha - a.alpha0) - d.clMax / a.CLa);
        const double qb = 0.5 * e.rho * V * V;
        const double L = qb * a.S * CL, D = qb * a.S * CD;

        // ---- autopilot
        double thCmd = 0, thr = 1.0;
        switch (ph) {
            case ROLL:
                thCmd = 0; thr = 1.0; intTh = 0;
                if (Vrx >= Vr) ph = ROTATE;
                break;
            case ROTATE:
                thCmd = rotTheta; thr = 1.0;
                if (z > 0.3) ph = CLIMB;
                break;
            case CLIMB:
                thr = 1.0;
                thCmd = clampd(0.12 + 0.06 * (V - Vclimb), 0.0, 0.25);
                if (z >= e.cruiseAlt - 4.0) { ph = CRUISE; tCruise = t; intAlt = 0; }
                break;
            case CRUISE: {
                const double eh = e.cruiseAlt - z;
                intAlt = clampd(intAlt + eh * dt, -4.0, 4.0);
                thCmd = clampd(0.05 + 0.04 * eh - 0.10 * vz + 0.02 * intAlt, -0.15, 0.20);
                thr = clampd(0.35 + 0.25 * (a.Vcruise - V), 0.05, 1.0);
                break;
            }
        }
        const double Vs = std::max(V, 0.6 * a.Vstall);
        const double B = 0.5 * e.rho * Vs * Vs * a.S * d.chord * std::fabs(a.Cmde) / a.Iyy;
        const double wc = 8.0;
        const double Kp = wc * wc / B, Kq = 2 * 0.8 * wc / B, Ki = 0.5 * Kp;
        const double err = thCmd - th;
        // positive elevator = trailing edge down = nose-DOWN moment, so pitch-up needs negative deflection
        const double deRaw = -(Kp * err - Kq * q + Ki * intTh);
        const double de = clampd(deRaw, -deMax, deMax);
        if (std::fabs(deRaw) < deMax) intTh = clampd(intTh + err * dt, -1.0, 1.0);

        // ---- forces
        const double T = thr * d.thrustStatic * std::max(0.0, 1.0 - V / d.vPitch);
        double Fx = T * std::cos(th) - D * std::cos(gam) - L * std::sin(gam);
        double Fz = T * std::sin(th) - D * std::sin(gam) + L * std::cos(gam) - d.mass * G;
        const double Cm = a.Cm0 + a.Cma * alpha + a.Cmq * q * d.chord / (2 * V) + a.Cmde * de;
        double M = qb * a.S * d.chord * Cm;

        const bool onGround = (z <= 0.0 && vz <= 0.0 && Fz <= 0.0);
        if (onGround) {
            const double N = -Fz;
            Fz = 0; vz = 0; z = 0;
            if (vx > 0.001) Fx -= e.rollFriction * N;
            else if (Fx <= e.rollFriction * N) { Fx = 0; vx = 0; }
            else Fx -= e.rollFriction * N;
            if (th <= 0.0 && M < 0.0) { M = 0; q = 0; th = 0; }       // nose on ground
            if (th >= d2r(15.0) && M > 0.0) { M = 0; q = 0; th = d2r(15.0); } // tail strike limit
        }

        const double n = L / (d.mass * G);

        // ---- integrate (semi-implicit Euler)
        vx += Fx / d.mass * dt;
        vz += Fz / d.mass * dt;
        x += vx * dt;
        z += vz * dt;
        q += M / a.Iyy * dt;
        th += q * dt;
        if (z < 0) { z = 0; if (vz < 0) vz = 0; }

        // ---- bookkeeping / failure checks
        if (!std::isfinite(x + z + th + q + vx + vz)) { R.reason = "numerical blow-up (unstable)"; break; }
        R.maxN = std::max(R.maxN, n);
        R.maxAlt = std::max(R.maxAlt, z);
        if (!R.lifted && z > 0.05) {
            R.lifted = true; tLift = t; R.rollDist = x; R.liftoffV = V; R.liftoffT = t;
            R.crabDeg = r2d(std::asin(std::min(1.0, std::fabs(cw) / std::max(V, 1.0))));
        }
        if (!R.lifted && x > e.runwayLen) { R.reason = "ran out of runway before liftoff"; break; }
        if (R.lifted && z > 0.3) {
            const double margin = a.alphaStall - alpha;
            R.minAlphaMargin = std::min(R.minAlphaMargin, margin);
            stallRun = (margin < 0) ? stallRun + dt : 0.0;
            if (stallRun > 0.3) { R.reason = "stalled in flight"; break; }
        }
        if (R.lifted && z <= 0.0 && t - tLift > 0.5) { R.reason = "ground strike after liftoff"; break; }
        if (R.maxN > e.nLimit) { R.reason = "load factor limit exceeded"; break; }
        if (ph == CRUISE) {
            R.reachedCruise = true;
            if (R.tToAlt == 0) R.tToAlt = t;
            if (t - tCruise > 6.0) { const double ee = e.cruiseAlt - z; sumErr2 += ee * ee; ++nErr; }
        }
        if (wantTrace && t >= nextTrace) {
            R.trace.push_back({t, x, z, V, r2d(alpha), r2d(th), n, thr, r2d(de)});
            nextTrace += 1.0;
        }
    }
    if (R.reason.empty()) {
        if (!R.lifted) R.reason = "never lifted off";
        else if (!R.reachedCruise) R.reason = "did not reach cruise altitude in time";
        else if (nErr > 0) { R.altRMS = std::sqrt(sumErr2 / nErr); R.ok = true; }
        else R.reason = "no settled cruise data";
    }
    return R;
}

// ------------------------------------------------------------------ evaluation
struct Scenario { std::string name; Env env; double weight; unsigned seed; };

std::vector<Scenario> makeScenarios(const Env& user) {
    Env calm = user; calm.windSpeed = 0;
    Env severe = user;
    severe.windSpeed = std::max(user.windSpeed * 1.5, 5.0);
    severe.gustFrac = std::max(user.gustFrac * 1.5, 0.25);
    return {{"calm", calm, 0.2, 11}, {"as entered", user, 0.5, 22}, {"severe (1.5x wind/gust)", severe, 0.3, 33}};
}

double scoreOne(const Design& d, const Env& e, const SimResult& r) {
    const Aero a = computeAero(d, e.rho);
    if (!r.ok) {
        double c = 1000.0 + 200.0 * (1.0 - std::min(1.0, r.maxAlt / e.cruiseAlt));
        if (!r.lifted) c += 100.0;
        return c;
    }
    double c = 0;
    c += 20.0 * r.rollDist / e.runwayLen;                       // short takeoff
    c += 10.0 * a.Vstall / 6.0;                                 // low stall speed
    c += 15.0 * std::min(r.altRMS, 3.0);                        // holds altitude in turbulence
    c += 10.0 * 12.0 / std::max(r.LDcruise, 1.0);               // efficient cruise
    c += 40.0 * std::max(0.0, 0.087 - r.minAlphaMargin);        // >=5 deg stall margin
    c += 10.0 * r.crabDeg / 30.0;                               // crosswind handling
    c += 5.0 * r.maxN / e.nLimit;
    c += 80.0 * std::max(0.0, 0.08 - a.SM) + 80.0 * std::max(0.0, a.SM - 0.22); // static margin
    c += 150.0 * std::max(0.0, 0.32 - a.Vh); // keep enough tail volume for pitch authority/damping (tail weight & drag are not modelled)
    return c;
}

double evaluate(const Design& d, const Env& user, bool verbose = false) {
    double total = 0;
    for (const auto& s : makeScenarios(user)) {
        SimResult r = simulate(d, s.env, s.seed);
        const double c = scoreOne(d, s.env, r);
        total += s.weight * c;
        if (verbose) {
            std::cout << "  [" << std::left << std::setw(26) << s.name << "] ";
            if (!r.ok) std::cout << "FAIL: " << r.reason << "\n";
            else std::cout << std::fixed << std::setprecision(2)
                           << "roll " << r.rollDist << " m | liftoff " << r.liftoffV << " m/s"
                           << " | alt-RMS " << r.altRMS << " m | max n " << r.maxN
                           << " | stall margin " << r.minAlphaMargin * 180 / PI << " deg"
                           << " | crab " << r.crabDeg << " deg\n";
        }
    }
    return total;
}

// ------------------------------------------------------------------ optimizer
struct OptCfg {
    int randomSamples = 600, refineIters = 1200;
    double spanPct = 20, chordPct = 25, tailAreaPct = 50, tailArmPct = 30, thrustPct = 0;
    double cgLo = 0.18, cgHi = 0.38;
};

struct Vars { double v[6]; }; // span, chord, tailArea, tailArm, cg, thrust

Design applyVars(Design d, const Vars& x) {
    d.span = x.v[0]; d.chord = x.v[1]; d.tailArea = x.v[2];
    d.tailArm = x.v[3]; d.cgFrac = x.v[4]; d.thrustStatic = x.v[5];
    return d;
}

Design optimize(const Design& base, const Env& env, const OptCfg& cfg) {
    double lo[6], hi[6];
    const double c[6] = {base.span, base.chord, base.tailArea, base.tailArm, base.cgFrac, base.thrustStatic};
    const double p[5] = {cfg.spanPct, cfg.chordPct, cfg.tailAreaPct, cfg.tailArmPct, cfg.thrustPct};
    for (int i = 0; i < 4; ++i) { lo[i] = c[i] * (1 - p[i] / 100); hi[i] = c[i] * (1 + p[i] / 100); }
    lo[4] = cfg.cgLo; hi[4] = cfg.cgHi;
    lo[5] = c[5] * (1 - cfg.thrustPct / 100); hi[5] = c[5] * (1 + cfg.thrustPct / 100);
    hi[3] = std::min(hi[3], 0.65 * base.length); // tail must stay on the airframe

    std::mt19937 rng(2024);
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> N(0, 1);

    Vars best; for (int i = 0; i < 6; ++i) best.v[i] = c[i];
    double bestCost = evaluate(applyVars(base, best), env);
    std::cout << "Baseline cost: " << std::fixed << std::setprecision(2) << bestCost << "\n";

    const int total = cfg.randomSamples + cfg.refineIters;
    auto progress = [&](int i) {
        if (total >= 10 && i % (total / 10) == 0)
            std::cout << "  " << std::setw(3) << (100 * i / total) << "%  best cost " << bestCost << std::endl;
    };
    for (int i = 0; i < cfg.randomSamples; ++i) {
        Vars t; for (int k = 0; k < 6; ++k) t.v[k] = lo[k] + (hi[k] - lo[k]) * U(rng);
        const double cost = evaluate(applyVars(base, t), env);
        if (cost < bestCost) { bestCost = cost; best = t; }
        progress(i);
    }
    double step = 0.15;
    for (int i = 0; i < cfg.refineIters; ++i) {
        Vars t = best;
        for (int k = 0; k < 6; ++k)
            if (hi[k] > lo[k] && U(rng) < 0.6)
                t.v[k] = clampd(t.v[k] + step * (hi[k] - lo[k]) * N(rng), lo[k], hi[k]);
        const double cost = evaluate(applyVars(base, t), env);
        if (cost < bestCost) { bestCost = cost; best = t; }
        step = std::max(0.01, step * 0.997);
        progress(cfg.randomSamples + i);
    }
    std::cout << "  100%  best cost " << bestCost << "\n";
    return applyVars(base, best);
}

// ------------------------------------------------------------------ UI helpers
double ask(const std::string& label, double cur) {
    while (true) {
        std::cout << "  " << label << " [" << cur << "]: ";
        std::string s;
        if (!std::getline(std::cin, s)) return cur;
        if (s.empty()) return cur;
        try { return std::stod(s); } catch (...) { std::cout << "  (not a number)\n"; }
    }
}

void printDesign(const Design& d, const Env& e) {
    const Aero a = computeAero(d, e.rho);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\n--- AIRCRAFT ---\n"
              << "mass " << d.mass << " kg (payload " << d.payload << " kg = " << d.payload / 0.45359 << " lb)\n"
              << "span " << d.span << " m | chord " << d.chord << " m | area " << a.S << " m^2 | AR " << a.AR << "\n"
              << "wing loading " << d.mass / a.S << " kg/m^2 | length " << d.length << " m\n"
              << "tail area " << d.tailArea << " m^2 | tail arm " << d.tailArm << " m | tail volume " << a.Vh << "\n"
              << "CG " << d.cgFrac * 100 << "% chord | static margin " << a.SM * 100 << "% chord"
              << (a.SM < 0.05 ? "  <-- too unstable" : (a.SM > 0.25 ? "  <-- very nose-heavy/stiff" : "")) << "\n"
              << "static thrust " << d.thrustStatic << " N (T/W " << d.thrustStatic / (d.mass * G) << ") | prop pitch speed " << d.vPitch << " m/s\n"
              << "CD0 " << d.cd0 << " | CLmax " << d.clMax << " | CLalpha " << a.CLa << "/rad\n"
              << "stall speed " << a.Vstall << " m/s | rotate " << 1.15 * a.Vstall << " | climb " << 1.4 * a.Vstall
              << " | cruise " << a.Vcruise << " m/s\n"
              << "--- ENVIRONMENT ---\n"
              << "wind " << e.windSpeed << " m/s @ " << e.windAngleDeg << " deg off runway heading | gust intensity "
              << e.gustFrac * 100 << "% | rho " << e.rho << "\n"
              << "runway " << e.runwayLen << " m | friction " << e.rollFriction << " | cruise alt " << e.cruiseAlt << " m\n";
}

void editDesign(Design& d) {
    std::cout << "Press Enter to keep a value.\n";
    d.mass = ask("total mass kg (1.5 = worst case incl. payload)", d.mass);
    d.payload = ask("payload kg (info only)", d.payload);
    d.span = ask("wingspan m", d.span);
    d.chord = ask("wing chord m", d.chord);
    d.length = ask("fuselage length m", d.length);
    d.tailArea = ask("horizontal tail area m^2", d.tailArea);
    d.tailArm = ask("tail arm (CG to tail) m", d.tailArm);
    d.cgFrac = ask("CG location (fraction of chord)", d.cgFrac);
    d.thrustStatic = ask("static thrust N", d.thrustStatic);
    d.vPitch = ask("prop pitch speed m/s", d.vPitch);
    d.cd0 = ask("CD0 zero-lift drag", d.cd0);
    d.clMax = ask("CLmax", d.clMax);
}

void editEnv(Env& e) {
    std::cout << "Press Enter to keep a value.\n";
    e.windSpeed = ask("wind speed m/s (at 10 m)", e.windSpeed);
    e.windAngleDeg = ask("wind angle deg (0 = headwind, 90 = crosswind, 180 = tailwind)", e.windAngleDeg);
    e.gustFrac = ask("gust intensity (0.3 = 30% of wind speed)", e.gustFrac);
    e.rho = ask("air density kg/m^3", e.rho);
    e.runwayLen = ask("runway length m", e.runwayLen);
    e.rollFriction = ask("rolling friction coeff", e.rollFriction);
    e.cruiseAlt = ask("cruise altitude m", e.cruiseAlt);
    e.nLimit = ask("load factor limit g", e.nLimit);
}

void editOpt(OptCfg& o) {
    std::cout << "Press Enter to keep a value.\n";
    o.randomSamples = (int)ask("random-search samples", o.randomSamples);
    o.refineIters = (int)ask("local-refine iterations", o.refineIters);
    o.spanPct = ask("span variation +/- %", o.spanPct);
    o.chordPct = ask("chord variation +/- %", o.chordPct);
    o.tailAreaPct = ask("tail area variation +/- %", o.tailAreaPct);
    o.tailArmPct = ask("tail arm variation +/- %", o.tailArmPct);
    o.thrustPct = ask("motor thrust variation +/- % (0 = fixed motor)", o.thrustPct);
    o.cgLo = ask("CG min (fraction of chord)", o.cgLo);
    o.cgHi = ask("CG max (fraction of chord)", o.cgHi);
}

void runSingle(const Design& d, const Env& e) {
    printDesign(d, e);
    std::cout << "\n=== Scenario results ===\n";
    evaluate(d, e, true);

    SimResult r = simulate(d, e, 22, true);
    std::cout << "\n=== Time history, scenario as entered (1 s steps) ===\n"
              << " t(s)  x(m)  alt(m)  Vair  alpha  theta   n    thr  elev\n";
    std::cout << std::fixed << std::setprecision(1);
    for (const auto& row : r.trace)
        std::cout << std::setw(5) << row.t << std::setw(7) << row.x << std::setw(7) << row.z << std::setw(7) << row.V
                  << std::setw(7) << row.alphaDeg << std::setw(7) << row.thetaDeg << std::setw(6) << row.n
                  << std::setw(6) << row.thr << std::setw(7) << row.de << "\n";

    int pass = 0; const int N = 30;
    for (int s = 0; s < N; ++s) if (simulate(d, e, 1000 + s).ok) ++pass;
    std::cout << "\nMonte-Carlo (" << N << " random gust seeds, wind as entered): " << pass << "/" << N << " flights OK\n";
}

} // namespace

int main() {
    Design d; Env e; OptCfg o;
    std::cout << "RC Cessna 152 (1:6) flight simulator & optimizer\n";
    while (true) {
        std::cout << "\n=============== MENU ===============\n"
                     "1) Show / edit aircraft parameters\n"
                     "2) Edit wind & environment\n"
                     "3) Run simulation of current design\n"
                     "4) OPTIMIZE design (uses current parameters as baseline)\n"
                     "5) Edit optimizer limits\n"
                     "6) Show current parameters\n"
                     "0) Quit\n> ";
        std::string s;
        if (!std::getline(std::cin, s)) break;
        if (s == "0") break;
        else if (s == "1") { editDesign(d); printDesign(d, e); }
        else if (s == "2") { editEnv(e); }
        else if (s == "3") { runSingle(d, e); }
        else if (s == "4") {
            auto t0 = std::chrono::steady_clock::now();
            Design best = optimize(d, e, o);
            double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "\n=== OPTIMIZED DESIGN (" << std::fixed << std::setprecision(1) << secs << " s) ===\n";
            printDesign(best, e);
            std::cout << "\n=== Results of optimized design ===\n";
            evaluate(best, e, true);
            std::cout << "\nChange vs baseline: span " << std::showpos << std::setprecision(1)
                      << 100 * (best.span / d.span - 1) << "%, chord " << 100 * (best.chord / d.chord - 1)
                      << "%, tail area " << 100 * (best.tailArea / d.tailArea - 1) << "%, tail arm "
                      << 100 * (best.tailArm / d.tailArm - 1) << "%" << std::noshowpos << "\n";
            std::cout << "Adopt optimized design as the new baseline? (y/N): ";
            std::string a; std::getline(std::cin, a);
            if (!a.empty() && (a[0] == 'y' || a[0] == 'Y')) d = best;
        }
        else if (s == "5") editOpt(o);
        else if (s == "6") printDesign(d, e);
        else std::cout << "Unknown option.\n";
    }
    return 0;
}