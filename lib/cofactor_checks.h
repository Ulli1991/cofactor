// cofactor: check rows, bars, exclusions and the report.
// Row format (unchanged from the 2026-11 audits): "<name padded to 92> n=<count> max_rel=<worst>  at <where>", then " | <verdict>".
// Bars are read at run time from a file of lines "<threshold> [@tag] <row-name prefix>" (longest matching prefix wins; '#'
// comments); a bar with @tag applies only when the config name contains tag, and then its row must exist;
// a row without a bar is reported, not judged. Rows named "EXCLUDED: ..." are counted and listed, never judged.
#ifndef COFACTOR_CHECKS_H
#define COFACTOR_CHECKS_H
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

struct Stat { double max = 0; std::string where; long n = 0; };
static std::map<std::string, Stat> g_stat;

static void rec(const std::string &q, double err, const std::string &where)
{
    Stat &s = g_stat[q]; s.n++;
    if(!(err <= s.max)) {s.max = err; s.where = where;} /* NaN is always recorded as the worst */
}
/* NaN anywhere makes the error NaN (fmax would drop it), and rec() records NaN as the worst value */
static double relerr(const double *a, const double *b, int n) /* max|a-b| / max|b| */
{
    double d = 0, s = 0;
    for(int q = 0; q < n; q++) {if(std::isnan(a[q]) || std::isnan(b[q])) {return NAN;} d = fmax(d, fabs(a[q] - b[q])); s = fmax(s, fabs(b[q]));}
    return (s > 0) ? d / s : d;
}
static double relerr_s(const double *a, const double *b, const double *s, int n) /* max|a-b| / max(|b|, s) */
{
    double d = 0, m = 0;
    for(int q = 0; q < n; q++) {if(std::isnan(a[q]) || std::isnan(b[q]) || std::isnan(s[q])) {return NAN;} d = fmax(d, fabs(a[q] - b[q])); m = fmax(m, fmax(fabs(b[q]), s[q]));}
    return (m > 0) ? d / m : d;
}
static void exclude(const std::string &reason, const std::string &where) {rec("EXCLUDED: " + reason, 1, where);}

/* bars: threshold per row-name prefix */
static std::vector<std::pair<std::string, double>> g_bars;
static bool load_bars(const char *path, const std::string &cfg)
{
    FILE *f = fopen(path, "r"); if(!f) {printf("cofactor: no bar file %s\n", path); return false;}
    char line[1024];
    while(fgets(line, sizeof(line), f))
    {
        if(line[0] == '#' || line[0] == '\n') {continue;}
        char *e; double t = strtod(line, &e); if(e == line) {continue;}
        while(*e == ' ' || *e == '\t') {e++;}
        if(*e == '@')
        {
            char *t = e + 1; while(*e && *e != ' ' && *e != '\t') {e++;}
            std::string tag(t, e - t); while(*e == ' ' || *e == '\t') {e++;}
            if(cfg.find(tag) == std::string::npos) {continue;}
        }
        std::string p(e); while(!p.empty() && (p.back() == '\n' || p.back() == '\r')) {p.pop_back();}
        if(!p.empty()) {g_bars.push_back({p, t});}
    }
    fclose(f); return true;
}
static const std::pair<std::string, double> *bar_for(const std::string &row)
{
    const std::pair<std::string, double> *b = nullptr;
    for(auto &x : g_bars) {if(row.compare(0, x.first.size(), x.first) == 0 && (!b || x.first.size() > b->first.size())) {b = &x;}}
    return b;
}

/* prints every row with its verdict and a summary; returns 0 if every barred row passes and every bar found its row,
   1 on any FAIL, 2 if nothing failed but the run is incomplete (forward-only build) */
static int cofactor_report(const char *target)
{
    int fail = 0, pass = 0, unbarred = 0; long excluded = 0;
#ifdef FWD_ONLY
    g_stat["REVERSE SKIPPED: forward-only build (reverse-mode rows absent)"].n = 1;
#endif
    for(auto &kv : g_stat)
    {
        const Stat &s = kv.second; const char *v; char buf[64];
        if(kv.first.compare(0, 9, "EXCLUDED:") == 0) {excluded += s.n; v = "excluded";}
        else if(kv.first.compare(0, 16, "REVERSE SKIPPED:") == 0) {v = "NOT RUN";}
        else {
            const std::pair<std::string, double> *b = bar_for(kv.first);
            if(!b) {unbarred++; v = "-";}
            else if(s.max <= b->second) {pass++; snprintf(buf, sizeof(buf), "PASS (bar %.0e)", b->second); v = buf;}
            else {fail++; snprintf(buf, sizeof(buf), "FAIL (bar %.0e)", b->second); v = buf;}
        }
        printf("%-92s n=%6ld max_rel=%.3e  at %s | %s\n", kv.first.c_str(), s.n, s.max, s.where.c_str(), v);
    }
    for(auto &x : g_bars)
    {
        bool hit = false; for(auto &kv : g_stat) {if(kv.first.compare(0, x.first.size(), x.first) == 0) {hit = true; break;}}
        if(!hit) {printf("cofactor: BAR WITHOUT ROW: %s (the check did not run)\n", x.first.c_str()); fail++;}
    }
    int rev = 1;
#ifdef FWD_ONLY
    rev = 0;
#endif
    printf("cofactor %s: %d PASS, %d FAIL, %d rows without bar, %ld records excluded, reverse mode %s\n",
           target, pass, fail, unbarred, excluded, rev ? "run" : "SKIPPED");
    return fail ? 1 : (rev ? 0 : 2);
}
#endif
