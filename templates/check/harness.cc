// cofactor check harness: calls the copied code (ext_*.h) and records checks; see docs/WRITING_A_CHECK.md.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "cofactor_enzyme.h"   // __enzyme_fwddiff, __enzyme_autodiff, enzyme_const, enzyme_dup
#include "cofactor_checks.h"   // rec(), relerr(), exclude(), load_bars(), cofactor_report()
#include "ext_code.h"          // written by spec.py

int main(int argc, char **argv)
{
    std::string cfg = (argc > 1) ? argv[1] : "default";
    if(!getenv("COFACTOR_BARS") || !load_bars(getenv("COFACTOR_BARS"), cfg)) {printf("cofactor: set COFACTOR_BARS\n"); return 4;}
    for(int i = 0; i < 1000; i++)
    {
        double x = 0.1 + i * 1e-3;
        char w[64]; snprintf(w, sizeof(w), "%s x=%g", cfg.c_str(), x);
        // row 0: the copied code reproduces what the original produces (here: just finite)
        double f = my_function(x);
        rec("0 copied code runs (1 = not finite)", std::isfinite(f) ? 0 : 1, w);
        // a derivative check: compare a hand-coded derivative, an identity or a booking with Enzyme
        double d = __enzyme_fwddiff<double>((void *)my_function, x, 1.0);
        rec("(1) TODO: hand derivative = Enzyme", fabs(d - d) / fmax(fabs(d), 1e-300), w);
#ifndef FWD_ONLY
        double r = __enzyme_autodiff<double>((void *)my_function, x);
        rec("(2) Enzyme reverse = forward", fabs(r - d) / fmax(fabs(d), 1e-300), w);
#endif
    }
    return cofactor_report("TODO_check_name");
}
