/* cooling.c: a small piece of "user code" for the cofactor example (any C/C++ code works the same way).
   A tabulated cooling function Lambda(T), log-log interpolated, its hand-coded derivative dLambda/dT, and the residual of a
   backward-Euler cooling step with the hand-coded derivative a Newton iteration would use. */
#include <math.h>

static const double T_edges[5] = {1.0e2, 1.0e4, 1.0e5, 1.0e6, 1.0e8};
static const double L_at_edge[5] = {1.0e-27, 1.0e-24, 5.0e-22, 1.0e-22, 3.0e-23};

/* Lambda(T) [erg cm^3 s^-1]: power law between tabulated points, constant outside the table */
double cooling_lambda(double T)
{
    if(T <= T_edges[0]) {return L_at_edge[0];}
    if(T >= T_edges[4]) {return L_at_edge[4];}
    int k = 0;
    while(T > T_edges[k + 1]) {k++;}
    double s = log(L_at_edge[k + 1] / L_at_edge[k]) / log(T_edges[k + 1] / T_edges[k]);
    return L_at_edge[k] * pow(T / T_edges[k], s);
}

/* dLambda/dT, hand-coded */
double cooling_dlambda_dT(double T)
{
    if(T <= T_edges[0] || T >= T_edges[4]) {return 0;}
    int k = 0;
    while(T > T_edges[k + 1]) {k++;}
    double s = log(L_at_edge[k + 1] / L_at_edge[k]) / log(T_edges[k + 1] / T_edges[k]);
#ifdef PLANTED_BUG
    if(k == 2) {s = log(L_at_edge[k] / L_at_edge[k - 1]) / log(T_edges[k] / T_edges[k - 1]);} /* demo: slope of the neighbouring segment */
#endif
    return s * cooling_lambda(T) / T;
}

/* backward-Euler cooling step in temperature units: R(T) = T - T0 + n2dt Lambda(T) = 0 */
double cooling_residual(double T, double T0, double n2dt)
{
    return T - T0 + n2dt * cooling_lambda(T);
}

/* dR/dT, hand-coded, as used by the Newton iteration */
double cooling_residual_dT(double T, double T0, double n2dt)
{
    (void)T0;
    return 1.0 + n2dt * cooling_dlambda_dT(T);
}
