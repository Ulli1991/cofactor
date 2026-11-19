// cofactor example (GIZMO): softened gravity kernel_gravity (mesh/kernel.h), copied verbatim. At h = 1, with
// phi(u) = mode -1, m1 = mode 1, m2 = mode 2, m0 = mode 0 (GIZMO's definitions, kernel.h comment):
//  (g1) m1 = (1/u) dphi/du          (g2) m2 = -(1/u) dm1/du          (g3) m0 = -(phi + u dphi/du)
//  (g4) Poisson: 3 m1 + u dm1/du = 4 pi W(u), W = kernel_main of the same KERNEL_FUNCTION (3D, normalised)
//  (g5) u -> 1 from below: phi = -1, m1 = 1, m2 = 3, m0 = 0 (Newtonian, as kernel_gravity returns for u >= 1)
//  (g6) continuity of phi, m1, m2, m0 across the internal joints of the piecewise kernels
// Errors are relative to the largest |reference| of the row over u in (0,1) (floor 1e-3 of it at zero crossings).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "cofactor_enzyme.h"
#include "cofactor_checks.h"
#include "ext_kernel_defs.h"
#include "ext_kernel_fun.h"

static double g(double u, int mode) {return kernel_gravity(u, 1.0, 1.0, mode);}
static double phi(double u) {return g(u, -1);}
static double m1(double u) {return g(u, 1);}
static double W(double u) {double wk = 0, dwk = 0; kernel_main(u, 1.0, 1.0, &wk, &dwk, -1); return wk;}

int main(int argc, char **argv)
{
    std::string cfg = (argc > 1) ? argv[1] : "default";
    if(!getenv("COFACTOR_BARS") || !load_bars(getenv("COFACTOR_BARS"), cfg)) {printf("cofactor: set COFACTOR_BARS\n"); return 4;}
    const int N = 4000;
    std::vector<double> us; for(int i = 0; i < N; i++) {us.push_back((i + 0.5) / N);}
    double s1 = 0, s2 = 0, s0 = 0, sp = 0;
    for(double u : us) {s1 = fmax(s1, fabs(m1(u))); s2 = fmax(s2, fabs(g(u, 2))); s0 = fmax(s0, fabs(g(u, 0))); sp = fmax(sp, fabs(4 * M_PI * W(u)));}
    for(double u : us)
    {
        char w[64]; snprintf(w, sizeof(w), "%s u=%.6f", cfg.c_str(), u);
        double dphi = __enzyme_fwddiff<double>((void *)phi, u, 1.0);
        double dm1 = __enzyme_fwddiff<double>((void *)m1, u, 1.0);
        double a1 = m1(u), b1 = dphi / u;
        rec("(g1) mode 1 = (1/u) dphi/du", fabs(a1 - b1) / fmax(fabs(b1), 1e-3 * s1), w);
        double a2 = g(u, 2), b2 = -dm1 / u;
        rec("(g2) mode 2 = -(1/u) d(mode 1)/du", fabs(a2 - b2) / fmax(fabs(b2), 1e-3 * s2), w);
        double a0 = g(u, 0), b0 = -(phi(u) + u * dphi);
        rec("(g3) mode 0 = -(phi + u dphi/du)", fabs(a0 - b0) / fmax(fabs(b0), 1e-3 * s0), w);
        double lhs = 3 * a1 + u * dm1, rhs = 4 * M_PI * W(u);
        rec("(g4) Poisson: 3 m1 + u dm1/du = 4 pi W (scale 4 pi max W)", fabs(lhs - rhs) / sp, w);
    }
    {
        double u = nextafter(1.0, 0.0);
        char w[64]; snprintf(w, sizeof(w), "%s u=1-", cfg.c_str());
        rec("(g5) Newtonian at u = 1: |phi + 1|", fabs(phi(u) + 1), w);
        rec("(g5) Newtonian at u = 1: |m1 - 1|", fabs(m1(u) - 1), w);
        rec("(g5) Newtonian at u = 1: |m2 - 3| / 3", fabs(g(u, 2) - 3) / 3, w);
        rec("(g5) Newtonian at u = 1: |m0| / max|m0|", fabs(g(u, 0)) / s0, w);
    }
    std::vector<double> joints;
#if (KERNEL_FUNCTION == 3)
    joints = {0.5};
#elif (KERNEL_FUNCTION == 4)
    joints = {0.2, 0.6};
#elif (KERNEL_FUNCTION == 5)
    joints = {1. / 3., 2. / 3.};
#elif (KERNEL_FUNCTION == 8)
    joints = {KERNEL_U0};
#endif
    for(double j : joints)
    {
        double lo = nextafter(j, 0.0), hi = j;
        const char *nm[4] = {"phi", "m1", "m2", "m0"}; int md[4] = {-1, 1, 2, 0}; double sc[4] = {fabs(phi(0.5)), s1, s2, s0};
        for(int q = 0; q < 4; q++)
        {
            char w[96]; snprintf(w, sizeof(w), "%s %s at u=%.6f: left %.17g right %.17g", cfg.c_str(), nm[q], j, g(lo, md[q]), g(hi, md[q]));
            rec(std::string("(g6) continuity at internal joints: ") + nm[q], fabs(g(lo, md[q]) - g(hi, md[q])) / sc[q], w);
        }
    }
    if(joints.empty()) {rec("(g6) continuity at internal joints: none in this kernel", 0, cfg);}
    return cofactor_report("gizmo/gravity_kernel");
}
