// cofactor example c_cooling: hand-coded derivatives of a tabulated cooling function, checked against Enzyme.
// The functions are the real ones from cooling.c (copied by spec.py); the harness only calls them.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "cofactor_enzyme.h"
#include "cofactor_checks.h"
#include "ext_table.h"
#include "ext_cooling.h"

int main(int argc, char **argv)
{
    std::string cfg = (argc > 1) ? argv[1] : "default";
    if(!getenv("COFACTOR_BARS") || !load_bars(getenv("COFACTOR_BARS"), cfg)) {printf("cofactor: set COFACTOR_BARS to the bar file\n"); return 4;}
    const int N = 4000;
    for(int i = 0; i < N; i++)
    {
        double T = pow(10., 1.5 + 7.8 * (i + 0.5) / N); /* 30 K .. 2e9 K, never exactly on a table edge (kinks) */
        char w[64]; snprintf(w, sizeof(w), "%s T=%.6e", cfg.c_str(), T);
        /* 0: the copied code runs and is finite */
        double L = cooling_lambda(T), dL = cooling_dlambda_dT(T);
        rec("0 copied code: Lambda(T), dLambda/dT finite (1 = not)", (std::isfinite(L) && std::isfinite(dL)) ? 0 : 1, w);
        /* 1: hand dLambda/dT vs Enzyme forward mode */
        double e = __enzyme_fwddiff<double>((void *)cooling_lambda, T, 1.0);
        rec("(1) hand dLambda/dT = Enzyme (relative)", fabs(dL - e) / fmax(fabs(e), 1e-300), w);
        /* 2: Newton derivative of the backward-Euler residual */
        double T0 = 2 * T, n2dt = 1e20;
        double hR = cooling_residual_dT(T, T0, n2dt);
        double eR = __enzyme_fwddiff<double>((void *)cooling_residual, T, 1.0, enzyme_const, T0, enzyme_const, n2dt);
        rec("(2) hand dR/dT (Newton) = Enzyme (relative)", fabs(hR - eR) / fmax(fabs(eR), 1e-300), w);
#ifndef FWD_ONLY
        /* 3: reverse mode agrees with forward mode */
        double r = __enzyme_autodiff<double>((void *)cooling_lambda, T);
        rec("(3) Enzyme reverse = forward on Lambda(T)", fabs(r - e) / fmax(fabs(e), 1e-300), w);
#endif
    }
    return cofactor_report("c_cooling");
}
