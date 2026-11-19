// Tests of lib/cofactor_checks.h (no Enzyme needed): bars, tags, NaN, missing rows, exclusions, exit codes.
// Build twice: g++ -std=c++14 -Ilib tests/test_checks.cc            (expects exit codes 0/1)
//              g++ -std=c++14 -Ilib -DFWD_ONLY tests/test_checks.cc (a passing run must give 2: incomplete)
#include <cmath>
#include <cstdio>
#include <string>
#include <unistd.h>
#include "cofactor_checks.h"

static int g_failures = 0;
static std::string bars_file(const char *content)
{
    std::string p = std::string("/tmp/cofactor_test_bars_") + std::to_string((long)getpid()) + ".txt";
    FILE *f = fopen(p.c_str(), "w"); fputs(content, f); fclose(f); return p;
}
static void reset() {g_stat.clear(); g_bars.clear();}
static void expect(const char *name, int got, int want)
{
    if(got != want) {printf("TEST FAIL: %s: got %d, want %d\n", name, got, want); g_failures++;}
    else {printf("TEST ok: %s\n", name);}
}
#ifdef FWD_ONLY
static const int PASS_CODE = 2;
#else
static const int PASS_CODE = 0;
#endif

int main()
{
    std::string b;

    reset(); b = bars_file("1e-12 (1)\n"); load_bars(b.c_str(), "cfg");
    rec("(1) check", 1e-14, "x"); expect("within bar", cofactor_report("t"), PASS_CODE);

    reset(); load_bars(b.c_str(), "cfg");
    rec("(1) check", 1e-6, "x"); expect("above bar", cofactor_report("t"), 1);

    reset(); load_bars(b.c_str(), "cfg");
    rec("(1) check", NAN, "x"); expect("NaN fails", cofactor_report("t"), 1);

    reset(); load_bars(b.c_str(), "cfg");
    rec("(9) other row", 1e-20, "x"); expect("bar without its row fails", cofactor_report("t"), 1);

    reset(); b = bars_file("1e-12 @psi (2) psi row\n"); load_bars(b.c_str(), "mfm_ded");
    expect("tagged bar ignored for other configs", cofactor_report("t"), PASS_CODE);

    reset(); load_bars(b.c_str(), "mfm_psi");
    expect("tagged bar required for matching config", cofactor_report("t"), 1);

    reset(); b = bars_file("1e-12 (1)\n"); load_bars(b.c_str(), "cfg");
    rec("(1) check", 0, "x"); rec("info: unbarred", 1e3, "x"); exclude("reached an abort path", "y");
    expect("unbarred rows and exclusions are not judged", cofactor_report("t"), PASS_CODE);

    reset(); b = bars_file("1e-12 (1)\n1e-6 (1) long\n"); load_bars(b.c_str(), "cfg");
    rec("(1) long prefix row", 1e-8, "x"); rec("(1) short", 1e-14, "x");
    expect("longest matching prefix wins", cofactor_report("t"), PASS_CODE);

    double a[2] = {1, NAN}, c[2] = {1, 2};
    expect("relerr propagates NaN", std::isnan(relerr(a, c, 2)) ? 1 : 0, 1);
    double d[2] = {1, 2.000002};
    expect("relerr value", std::fabs(relerr(d, c, 2) - 1e-6) < 1e-12 ? 1 : 0, 1);

    printf("%s: %d failure(s)\n", g_failures ? "FAILED" : "PASSED", g_failures);
    return g_failures ? 1 : 0;
}
