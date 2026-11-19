// cofactor target example_kernel: GIZMO's kernel_main (mesh/kernel.h), copied verbatim.
// Checks, for one KERNEL_FUNCTION and NUMDIMS (compile-time defines of the config):
//   row 0: mode 0 (wk and dwk together) gives bitwise the values of mode -1 (wk only) and mode +1 (dwk only);
//   (k)  : the hand-coded dwk equals Enzyme's d(wk)/du times hinv4/hinv3 (dwk carries hinv4, wk carries hinv3);
//   (k)  : Enzyme reverse mode equals forward mode.
// This is the smallest complete target: use it as the template for new ones and as the smoke test of an installation.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "cofactor_enzyme.h"
#include "cofactor_checks.h"
#include "ext_kernel_defs.h"
#include "ext_kernel_fun.h"

/* the function Enzyme differentiates: wk(u) at fixed h */
static double wk_of_u(double u, double hinv3, double hinv4)
{
    double wk = 0, dwk = 0;
    kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
    return wk;
}

int main(int argc, char **argv)
{
    std::string cfg = (argc > 1) ? argv[1] : "default";
    if(!getenv("COFACTOR_BARS") || !load_bars(getenv("COFACTOR_BARS"), cfg)) {printf("cofactor: set COFACTOR_BARS to the bar file\n"); return 4;}
    printf("example_kernel %s: KERNEL_FUNCTION=%d NUMDIMS=%d KERNEL_NORM=%.17g\n", cfg.c_str(), (int)KERNEL_FUNCTION, (int)NUMDIMS, (double)KERNEL_NORM);
    const double hs[3] = {1.0, 0.37, 2.9};
    const int N = 2000;
    for(double h : hs)
    {
        double hinv, hinv3, hinv4; kernel_hinv(h, &hinv, &hinv3, &hinv4);
        double scale = KERNEL_NORM * hinv4; /* dwk = (dw/du) KERNEL_NORM hinv4, dw/du of order 1-20 */
        for(int i = 0; i <= N + 2; i++)
        {
            double u = (i < N) ? (i + 0.5) / N : ((i == N) ? 0.5 : 0.0); /* interior samples, the spline joint, the centre */
            char w[96]; snprintf(w, sizeof(w), "%s h=%g u=%.6f", cfg.c_str(), h, u);
            double w0 = 0, d0 = 0, w1 = 0, d1 = 0, dd = 0;
            kernel_main(u, hinv3, hinv4, &w0, &d0, 0);
            kernel_main(u, hinv3, hinv4, &w1, &dd, -1);
            kernel_main(u, hinv3, hinv4, &dd, &d1, +1);
            rec("0 kernel_main: mode 0 = modes -1/+1, bitwise (1 = differs)", (w0 == w1 && d0 == d1) ? 0 : 1, w);
            double fw = __enzyme_fwddiff<double>((void *)wk_of_u, u, 1.0, enzyme_const, hinv3, enzyme_const, hinv4);
            rec("(k) hand dwk = Enzyme d(wk)/du x hinv4/hinv3 (scale KERNEL_NORM hinv4)", fabs(d1 - fw * (hinv4 / hinv3)) / scale, w);
#ifndef FWD_ONLY
            double rv = __enzyme_autodiff<double>((void *)wk_of_u, u, enzyme_const, hinv3, enzyme_const, hinv4);
            rec("(k) Enzyme reverse = forward on wk(u)", fabs(rv - fw) / (KERNEL_NORM * hinv3), w); /* d(wk)/du carries KERNEL_NORM hinv3 */
#endif
        }
    }
    return cofactor_report("example_kernel");
}
