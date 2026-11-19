// cofactor target mhd_face (port of enzyme/mhd_audit/harness/face_harness.cc, 2026-11-19).
// Enzyme face harness for GIZMO MHD (Dedner GLM and psi/c_h cleaning). The face code is the real code, copied verbatim
// from the build copy by gen.py (see EXTRACTED.txt); P/CellP are macros for arrays passed in as arguments (Enzyme cannot
// shadow pointer globals). Inputs: face states dumped by TEST_MHD_DUMP_FACE_STATES from a 1-rank run of the same build.
#include <stdio.h>
#include <cstdio> /* before the fflush stub: <cstdio> #undef's function macros */
#include <string>
#include <map>
#include <vector>
#include <mpi.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../declarations/allvars.h"
#include "cofactor_glue.h" /* GIZMO warning/abort stubs (counted), Enzyme declarations, inactive globals/functions */
#include "ext_head.h"
__attribute__((used)) void *__enzyme_inactive_global_TBA = (void *)&TimeBinActive;

#define P Pw
#define CellP Cw
#define Get_Gas_CleaningSpeed(i) H_Get_Gas_CleaningSpeed(Pw, Cw, i)
#define Get_Gas_Psi_Conservative(i) H_Get_Gas_Psi_Conservative(Pw, Cw, i)
#define Get_Gas_PhiField(i) H_Get_Gas_PhiField(Pw, Cw, i)
#define estep_cell_mag_energy(i) H_estep_cell_mag_energy(Pw, Cw, i)
#include "ext_fun.h"

typedef struct INPUT_STRUCT_NAME LIn;
typedef struct OUTPUT_STRUCT_NAME LOut;
typedef struct particle_data PD;
typedef struct gas_cell_data CD;
struct FaceAux { double dt_j, fci, fcj; int jact, tbact; };
enum { O_D = 0, O_DC, O_WL, O_WR, O_JUMP, O_CORR, O_DBS, O_BR = 7, O_BL = 10, O_A = 13, O_AN = 16, O_DTBI = 17, O_ADJI = 20, O_DTBJ = 29,
       O_ADJJ = 32, O_DIVDISS = 41, O_DTBDISS = 42, O_HEAT = 45, O_HEATJ = 46, O_DTB = 47, O_ACC = 50, O_DTU = 53, O_DIVB = 54,
       O_PHIM = 55, O_PHIDB = 56, O_PHICORR = 57, O_PM = 60, O_SM = 61, O_BNR = 62, O_BNL = 63, O_DTPSIQ = 65, O_ADVHEAT = 66,
       O_MDOT = 67, O_REACHED = 68, O_LEDG1 = 69, O_FB = 70, O_FP = 73, O_FV = 74, O_PHIR = 77, O_PHIL = 78, O_JAC = 80, O_DTMASS = 104, O_DTPHI = 105, NOUT = 106 };

/* outputs, read element by element (no struct copies on the differentiated path) */
static void pack(const LOut &out, const struct Conserved_var_Riemann &F, const struct Riemann_outputs &Ro, const struct Input_vec_Riemann &Rv,
                 const double *A, double An, int reached, const double *jac, const CD *cj, double *v)
{
    for(int q = 0; q < NOUT; q++) {v[q] = 0;}
    v[O_D] = F.B_normal_corrected; v[O_DC] = v[O_D];
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    v[O_DC] = -An * Ro.psi_B_centered; v[O_WL] = Ro.psi_wL; v[O_WR] = Ro.psi_wR; v[O_JUMP] = Ro.psi_B_jump;
    v[O_CORR] = Ro.psi_corr; v[O_DBS] = Ro.psi_db_speed;
    for(int q = 0; q < 3; q++) {
        v[O_DTBI + q] = out.DtB_PsiCons[q]; v[O_DTBJ + q] = cj->DtB_PsiCons[q];
        for(int d = 0; d < 3; d++) {v[O_ADJI + 3 * q + d] = out.PsiGradAdj[q][d]; v[O_ADJJ + 3 * q + d] = cj->PsiGradAdj[q][d];}
        v[O_DTBDISS + q] = out.DtB_PsiDiss[q];
    }
    v[O_DIVDISS] = out.PsiDivDiss; v[O_HEAT] = out.PsiDissHeat; v[O_HEATJ] = cj->PsiDissHeat;
    v[O_DTPSIQ] = out.DtPsiQ; v[O_ADVHEAT] = out.PsiAdvHeat;
#endif
    for(int q = 0; q < 3; q++) {
        v[O_BR + q] = Rv.R.B[q]; v[O_BL + q] = Rv.L.B[q]; v[O_A + q] = A[q];
        v[O_DTB + q] = out.DtB[q]; v[O_ACC + q] = out.Acc[q]; v[O_FB + q] = F.B[q]; v[O_FV + q] = F.v[q];
#ifdef DIVBCLEANING_DEDNER_EVOLUTION
        v[O_PHICORR + q] = out.DtB_PhiCorr[q];
#endif
    }
    v[O_AN] = An; v[O_DTU] = out.DtInternalEnergy; v[O_DIVB] = out.divB;
    v[O_PHIM] = Ro.phi_normal_mean; v[O_PHIDB] = Ro.phi_normal_db; v[O_PM] = Ro.P_M; v[O_SM] = Ro.S_M;
    double bnr = 0, bnl = 0; if(An > 0) {for(int q = 0; q < 3; q++) {bnr += A[q] / An * Rv.R.B[q]; bnl += A[q] / An * Rv.L.B[q];}}
    v[O_BNR] = bnr; v[O_BNL] = bnl; v[O_MDOT] = F.rho; v[O_REACHED] = reached; v[O_FP] = F.p;
#ifdef TEST_MHD_ENERGY_LEDGER_OUTPUT
    v[O_LEDG1] = out.DtE_MHDLedger[1];
#endif
    v[O_PHIR] = Rv.R.phi; v[O_PHIL] = Rv.L.phi;
    for(int q = 0; q < 24; q++) {v[O_JAC + q] = jac[q];}
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    v[O_DTMASS] = out.DtMass;
#if defined(DIVBCLEANING_DEDNER_EVOLUTION)
    v[O_DTPHI] = out.DtPhi;
#endif
#endif
}

/* one evaluation of pair (i = local, j = index 1) exactly as hydro_force_evaluate does it; outputs packed into v */
static void face_eval(const LIn *Lin, PD *Pw, CD *Cw, const FaceAux *aux, double *vout)
{
#include "ext_decl.h"
    local = *Lin;
#include "ext_prelude.h"
    j = 1;
    dt_hydrostep_j = aux->dt_j; dt_hydrostep = DMAX(dt_hydrostep_i, dt_hydrostep_j);
    double FluxCorrectionFactor_to_i = aux->fci, FluxCorrectionFactor_to_j = aux->fcj;
    int j_is_active_for_fluxes = aux->jact;
    int reached = 0;
    do {
#include "ext_body.h"
        reached = 1;
    } while(0);
    double A3[3] = {Face_Area_Vec[0], Face_Area_Vec[1], Face_Area_Vec[2]}, jac[24];
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    for(int q = 0; q < 24; q++) {jac[q] = psi_face_jac[q / 8][q % 8];}
#else
    for(int q = 0; q < 24; q++) {jac[q] = 0;}
#endif
    pack(out, Fluxes, Riemann_out, Riemann_vec, A3, Face_Area_Norm, reached, jac, &Cw[1], vout);
    (void)k_hll; (void)b_hll; (void)v_hll;
}

/* final operations of hydro_toplevel.cc on cell index 1 */
static void fin_powell_convert(PD *Pw, CD *Cw)
{
    int i = 1;
#include "ext_fin_powell.h"
    double magnorm_closure = 1;
#include "ext_fin_convert.h"
    (void)magnorm_closure;
}
static void fin_psi(PD *Pw, CD *Cw)
{
    int i = 1;
#include "ext_fin_psi.h"
    (void)i;
}

/* the cell accumulators the pair loop and the final operations add to */
static void zero_acc(CD *c)
{
    c->divB = 0; c->DtB = {}; c->HydroAccel = {}; c->DtInternalEnergy = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    c->DtMass = 0;
#endif
#ifdef TEST_MHD_ENERGY_LEDGER_OUTPUT
    for(int q = 0; q < 5; q++) {c->DtE_MHDLedger[q] = 0;}
#endif
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    c->DtB_PsiCons = {}; c->DtB_PsiDiss = {}; for(int q = 0; q < 3; q++) {c->PsiGradAdj[q] = {};}
    c->PsiDivDiss = c->PsiDissHeat = c->DtPsiQ = c->PsiAdvHeat = 0; c->DtPsiZ = 0;
#endif
#ifdef DIVBCLEANING_DEDNER_EVOLUTION
    c->DtB_PhiCorr = {};
#endif
}

/* differentiable wrappers: the face writes into CellP[j], so they work on copies */
static void fvec(const LIn *Lin, PD *Pw, CD *Cw, const FaceAux *aux, double *v)
{
    PD P2[2]; CD C2[2]; memcpy(P2, Pw, sizeof(P2)); memcpy(C2, Cw, sizeof(C2));
    face_eval(Lin, P2, C2, aux, v);
}
/* the real reconstruction (reimann.h), alone: out2 = (Q_L, Q_R) */
static void recon_states(double Qi, double *gi, double Qj, double *gj, const double *di, const double *dj, int mode, double *out2)
{
    double d1[3] = {di[0], di[1], di[2]}, d2[3] = {dj[0], dj[1], dj[2]}, QL, QR;
    reconstruct_face_states(Qi, gi, Qj, gj, d1, d2, &QL, &QR, mode);
    out2[0] = QL; out2[1] = QR;
}
static double recon_comb(double *q, double *gi, double *gj, const double *di, const double *dj, int mode, double a, double b) /* a Q_L + b Q_R */
{
    double o[2]; recon_states(q[0], gi, q[1], gj, di, dj, mode, o); return a * o[0] + b * o[1];
}
static double fsel(const LIn *Lin, PD *Pw, CD *Cw, const FaceAux *aux, int which)
{
    double v[NOUT]; fvec(Lin, Pw, Cw, aux, v); return v[which];
}
static double Emag(PD *Pw, CD *Cw) {return H_estep_cell_mag_energy(Pw, Cw, 1);}
static double Etot(PD *Pw, CD *Cw) {return H_estep_cell_energy(Pw, Cw, 1);}

/* ---------------- reporting (format of the KETJU/RUMI audits) ---------------- */
#include "cofactor_checks.h" /* rec, relerr, relerr_s, exclude, bars, report */
static std::string g_cfg;
static double g_pw[8] = {0, 0, 0, 0, 0, 0, 0, 0}; /* Powell energy split sums */

/* ---------------- dump reader ---------------- */
struct Rec {
    LIn L; PD Pj; CD Cj0; LOut out0; double aux[4]; int ia[4];
    LOut out1; CD Cj1; struct Conserved_var_Riemann F; struct Riemann_outputs Ro; struct Input_vec_Riemann Rv; double jac[24]; double a4[4]; struct kernel_hydra kh;
};
static bool read_rec(FILE *f, Rec *r)
{
    size_t SA = sizeof(LIn) + sizeof(PD) + sizeof(CD) + sizeof(LOut) + 64;
    std::vector<char> a(SA); if(fread(a.data(), SA, 1, f) != 1) {return false;}
    char *b = a.data();
    memcpy(&r->L, b, sizeof(LIn)); b += sizeof(LIn); memcpy(&r->Pj, b, sizeof(PD)); b += sizeof(PD);
    memcpy(&r->Cj0, b, sizeof(CD)); b += sizeof(CD); memcpy(&r->out0, b, sizeof(LOut)); b += sizeof(LOut);
    memcpy(r->aux, b, 32); b += 32; memcpy(r->ia, b, 16);
    bool ok = fread(&r->out1, sizeof(LOut), 1, f) == 1 && fread(&r->Cj1, sizeof(CD), 1, f) == 1 && fread(&r->F, sizeof(r->F), 1, f) == 1 &&
              fread(&r->Ro, sizeof(r->Ro), 1, f) == 1 && fread(&r->Rv, sizeof(r->Rv), 1, f) == 1 && fread(r->jac, 8, 24, f) == 24 &&
              fread(r->a4, 8, 4, f) == 4 && fread(&r->kh, sizeof(r->kh), 1, f) == 1;
    return ok;
}

/* GIZMO's own values for this face, in pack() layout: accumulators as after - before */
/* vs: the scale each value can be recovered to; for the running sums GIZMO's per-face part is (after - before), known
   only to roundoff of the larger accumulated sum */
static void pack_dump(const Rec *r, double *v, double *vs)
{
    double v1[NOUT], v0[NOUT];
    pack(r->out1, r->F, r->Ro, r->Rv, r->a4, r->a4[3], 1, r->jac, &r->Cj1, v1);
    pack(r->out0, r->F, r->Ro, r->Rv, r->a4, r->a4[3], 1, r->jac, &r->Cj0, v0);
    for(int q = 0; q < NOUT; q++) {v[q] = v1[q]; vs[q] = fabs(v1[q]);}
    std::vector<int> acc = {O_DTBI, O_DTBI + 1, O_DTBI + 2, O_DTBJ, O_DTBJ + 1, O_DTBJ + 2, O_DIVDISS, O_DTBDISS, O_DTBDISS + 1, O_DTBDISS + 2, O_HEAT, O_HEATJ,
                 O_DTB, O_DTB + 1, O_DTB + 2, O_ACC, O_ACC + 1, O_ACC + 2, O_DTU, O_DIVB, O_PHICORR, O_PHICORR + 1, O_PHICORR + 2, O_DTPSIQ, O_ADVHEAT, O_LEDG1};
    for(int q = 0; q < 9; q++) {acc.push_back(O_ADJI + q); acc.push_back(O_ADJJ + q);}
    for(int q : acc) {v[q] = v1[q] - v0[q]; vs[q] = DMAX(fabs(v1[q]), fabs(v0[q]));}
}

/* derivative directions: 0-2 B_i, 3-5 B_j (intensive), 6-14 gradB_i[k][d], 15-23 gradB_j[k][d] */
static void seed(int dir, const CD *Cw1, LIn *dL, CD *dC1)
{
    if(dir < 3) {dL->BPred[dir] = 1;}
    else if(dir < 6) {dC1->BPred[dir - 3] = Cw1->Mass / Cw1->Density;}
    else if(dir < 15) {int q = dir - 6; dL->Gradients.B[q / 3][q % 3] = 1;}
    else {int q = dir - 15; dC1->Gradients.B[q / 3][q % 3] = 1;}
}
/* reverse-mode gradient in the same 24 components */
static void gather(const CD *Cw1, const LIn *bL, const CD *bC1, double *g)
{
    for(int k = 0; k < 3; k++) {g[k] = bL->BPred[k]; g[3 + k] = bC1->BPred[k] * Cw1->Mass / Cw1->Density;}
    for(int q = 0; q < 9; q++) {g[6 + q] = bL->Gradients.B[q / 3][q % 3]; g[15 + q] = bC1->Gradients.B[q / 3][q % 3];}
}

/* a cell for side i built from the hydro input of i, on the template of cell j (used for the per-cell energy checks) */
static void cell_from_input(const LIn *L, const PD *Pt, const CD *Ct, PD *Pi, CD *Ci)
{
    *Pi = *Pt; *Ci = *Ct;
    Pi->Mass = L->Mass; Ci->Mass = L->Mass; Ci->Density = L->Density; Pi->Type = 0;
    double V = L->Mass / L->Density;
    for(int q = 0; q < 3; q++) {Ci->BPred[q] = L->BPred[q] * V; Ci->B[q] = Ci->BPred[q]; Ci->VelPred[q] = L->Vel[q]; Pi->Vel[q] = L->Vel[q];}
    Ci->InternalEnergy = Ci->InternalEnergyPred = L->InternalEnergyPred; Ci->Pressure = L->Pressure;
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    Ci->PsiCleaningSpeed = L->PsiCh;
    Ci->PsiZPred = Ci->PsiZ = L->PsiCons * sqrt(All.cf_atime * V / All.cf_a3inv) / L->PsiCh;
#endif
}
static void sync_cell(PD *Pc, CD *Cc) /* stored = predicted state, as at a full sync */
{
    Cc->B = Cc->BPred; Cc->InternalEnergy = Cc->InternalEnergyPred; for(int q = 0; q < 3; q++) {Pc->Vel[q] = Cc->VelPred[q];}
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    Cc->PsiZ = Cc->PsiZPred;
#endif
    Pc->Type = 0;
}

/* final rates of cell i from this one face, through the real final operations (Powell, total->thermal conversion, psi/c_h
   block); R: 0 du/dt (incl. psi heat), 1-3 dv/dt, 4-6 d(BV)/dt, 7 dZ/dt, 8 ledger row 0 (Powell energy), 9 m du/dt */
enum { R_DU = 0, R_DV = 1, R_DB = 4, R_DZ = 7, R_LED0 = 8, R_MDU = 9, NR = 10 };
static void rates_i(const LIn *Lin, PD *Pw, CD *Cw, const FaceAux *aux, double *R)
{
    double v[NOUT]; PD P2[2]; CD C2[2]; memcpy(P2, Pw, sizeof(P2)); memcpy(C2, Cw, sizeof(C2));
    face_eval(Lin, P2, C2, aux, v);
    PD Pt[2]; CD Ct[2]; cell_from_input(Lin, &Pw[1], &Cw[1], &Pt[1], &Ct[1]); sync_cell(&Pt[1], &Ct[1]); zero_acc(&Ct[1]); Pt[0] = Pt[1]; Ct[0] = Ct[1];
    for(int k = 0; k < 3; k++) {Ct[1].HydroAccel[k] = v[O_ACC + k]; Ct[1].DtB[k] = v[O_DTB + k];}
    Ct[1].DtInternalEnergy = v[O_DTU]; Ct[1].divB = v[O_DIVB];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    Ct[1].DtMass = v[O_DTMASS];
#endif
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    for(int k = 0; k < 3; k++) {Ct[1].DtB_PsiCons[k] = v[O_DTBI + k]; Ct[1].DtB_PsiDiss[k] = v[O_DTBDISS + k];}
    Ct[1].PsiDivDiss = v[O_DIVDISS]; Ct[1].PsiDissHeat = v[O_HEAT]; Ct[1].DtPsiQ = v[O_DTPSIQ];
#endif
    fin_powell_convert(Pt, Ct);
    double heat = 0;
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
    fin_psi(Pt, Ct); heat = Ct[1].PsiDissHeat; R[R_DZ] = Ct[1].DtPsiZ;
#else
    R[R_DZ] = 0;
#endif
    R[R_MDU] = Ct[1].DtInternalEnergy + heat; R[R_DU] = R[R_MDU] / Pt[1].Mass; /* hydro_toplevel.cc: /= Mass, then += PsiDissHeat/Mass */
    for(int k = 0; k < 3; k++) {R[R_DV + k] = Ct[1].HydroAccel[k]; R[R_DB + k] = Ct[1].DtB[k];}
#ifdef TEST_MHD_ENERGY_LEDGER_OUTPUT
    R[R_LED0] = Ct[1].DtE_MHDLedger[0];
#else
    R[R_LED0] = 0;
#endif
}
/* Powell pieces of this face for cell i and the candidate heat bookings:
   T: 0-2 Powell momentum rate M_P (= m dv/dt), 3-5 Powell induction d(BV)/dt, 6 Powell energy (ledger row 0), 7 face divergence D,
   8 heat_a = -(v.M_P + B.dB_P) (Powell work into heat), 9 heat_c = same with v -> v - w_f (w_f = (v_i+v_j)/2, the face frame),
   10 -E_P + v.M_P (Janhunen u minus code u), 11 v.M_P + B.dB_P (Enzyme-checked kinetic+magnetic Powell work) */
enum { T_MP = 0, T_BP = 3, T_EP = 6, T_D = 7, T_HA = 8, T_HC = 9, T_JAN = 10, T_W = 11, NT = 12 };
static void powell_terms(const LIn *Lin, PD *Pw, CD *Cw, const FaceAux *aux, double *T)
{
    double v[NOUT]; PD P2[2]; CD C2[2]; memcpy(P2, Pw, sizeof(P2)); memcpy(C2, Cw, sizeof(C2));
    face_eval(Lin, P2, C2, aux, v);
    PD Pt[2]; CD Ct[2]; cell_from_input(Lin, &Pw[1], &Cw[1], &Pt[1], &Ct[1]); sync_cell(&Pt[1], &Ct[1]); zero_acc(&Ct[1]); Pt[0] = Pt[1]; Ct[0] = Ct[1];
    Ct[1].divB = v[O_D];
    fin_powell_convert(Pt, Ct); /* Powell only: HydroAccel = M_P/m, DtB = dB_P, DtInternalEnergy = m du of Powell */
    Vec3<double> Bi = Ct[1].Bfield() * All.cf_a2inv, vi = Ct[1].VelPred / All.cf_atime, vj = Cw[1].VelPred / All.cf_atime, wf = 0.5 * (vi + vj);
    double vM = 0, BdB = 0, rM = 0, rB = 0;
    for(int k = 0; k < 3; k++) {double M = Ct[1].HydroAccel[k] * Pt[1].Mass; T[T_MP + k] = M; T[T_BP + k] = Ct[1].DtB[k];
        vM += vi[k] * M; BdB += Bi[k] * Ct[1].DtB[k]; rM += (vi[k] - wf[k]) * M; rB += Bi[k] * Ct[1].DtB[k] * 0; }
    /* relative-velocity version of the induction work: dB_P = -D v_i, so B.dB_P with v_i -> v_i - w_f */
    double rBdB = 0; for(int k = 0; k < 3; k++) {rBdB += Bi[k] * (-v[O_D]) * (vi[k] - wf[k]);}
#ifdef TEST_MHD_ENERGY_LEDGER_OUTPUT
    T[T_EP] = Ct[1].DtE_MHDLedger[0];
#else
    T[T_EP] = 0;
#endif
    T[T_D] = v[O_D]; T[T_HA] = -(vM + BdB); T[T_HC] = -(rM + rBdB); T[T_JAN] = -T[T_EP] + vM; T[T_W] = vM + BdB; (void)rB;
}

/* uniform boost of every velocity the pair reads (cell i hydro input, j's VelPred, ParticleVel in MFV, P.Vel) along axis k */
static void boost_seed(int k, LIn *dL, PD *dP, CD *dC)
{
    memset(dL, 0, sizeof(LIn)); memset(dP, 0, 2 * sizeof(PD)); memset(dC, 0, 2 * sizeof(CD));
    dL->Vel[k] = 1; dC[1].VelPred[k] = 1; dP[1].Vel[k] = 1;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    dL->ParticleVel[k] = 1; dC[1].ParticleVel[k] = 1;
#endif
}
static double gmax(const double *a, int n) { double m = 0; for(int q = 0; q < n; q++) {m = DMAX(m, fabs(a[q]));} return m; }

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    if(argc < 3) {printf("usage: face_harness <face_dump.bin> <cfg> [show N]\n"); return 1;}
    g_cfg = argv[2];
    if(!getenv("COFACTOR_BARS") || !load_bars(getenv("COFACTOR_BARS"), g_cfg)) {printf("cofactor: set COFACTOR_BARS to the bar file\n"); return 1;}
    int show = (argc > 4) ? atoi(argv[4]) : -1;
    FILE *f = fopen(argv[1], "rb"); if(!f) {printf("no dump\n"); return 1;}
    int sz[10]; if(fread(sz, sizeof(sz), 1, f) != 1) {return 1;}
    int mine[10] = {0x475a4644, (int)sizeof(LIn), (int)sizeof(PD), (int)sizeof(CD), (int)sizeof(LOut), (int)sizeof(struct Conserved_var_Riemann),
                    (int)sizeof(struct Riemann_outputs), (int)sizeof(struct Input_vec_Riemann), (int)sizeof(struct kernel_hydra), (int)sizeof(All)};
    for(int q = 0; q < 10; q++) {if(sz[q] != mine[q]) {printf("layout mismatch at %d: dump %d harness %d\n", q, sz[q], mine[q]); return 2;}}
    if(fread(&All, sizeof(All), 1, f) != 1) {return 1;}
    double g4[4]; if(fread(g4, sizeof(g4), 1, f) != 1) {return 1;}
    fac_mu = g4[0]; fac_magnetic_pressure = g4[1]; boxSize = g4[2]; boxHalf = g4[3];
    printf("layout ok; cf_atime=%g a2inv=%g a3inv=%g comoving=%d\n", All.cf_atime, All.cf_a2inv, All.cf_a3inv, All.ComovingIntegrationOn);
    static Rec r; long nf = 0;
    long cat_count[8] = {0};
    while(read_rec(f, &r))
    {
        long id = nf++;
        char w[256]; snprintf(w, sizeof(w), "%s face=%ld step=%d", g_cfg.c_str(), id, r.ia[2]);
        PD Pw[2] = {r.Pj, r.Pj}; CD Cw[2] = {r.Cj0, r.Cj0};
        zero_acc(&Cw[1]); zero_acc(&Cw[0]);
        TimeBinActive[r.Pj.TimeBin] = r.ia[1];
        FaceAux aux = {r.aux[0], r.aux[1], r.aux[2], r.ia[0], r.ia[1]};
        double v[NOUT], vd[NOUT], vs[NOUT];
        g_warn_hits = g_abort_hits = 0;
        fvec(&r.L, Pw, Cw, &aux, v);
        pack_dump(&r, vd, vs);
        if(g_abort_hits) {rec("EXCLUDED: faces that reach a GIZMO abort (endrun) path", 1, w); continue;}
        if(g_warn_hits) {rec("EXCLUDED: faces that reach a GIZMO warning path", 1, w); continue;}
        if(v[O_REACHED] == 0) {rec("faces skipped by the pair filter (r2/h checks)", 0, w); continue;}
        /* 0: the harness reproduces GIZMO's own face */
        {
            double e1 = relerr(v + O_D, vd + O_D, 1), e2 = relerr(v + O_FB, vd + O_FB, 3), e3 = relerr_s(v + O_DTB, vd + O_DTB, vs + O_DTB, 3), e4 = relerr_s(v + O_ACC, vd + O_ACC, vs + O_ACC, 3);
            double e5 = relerr_s(v + O_DTU, vd + O_DTU, vs + O_DTU, 1), e6 = relerr(v + O_BR, vd + O_BR, 6);
            rec("0 harness = GIZMO run: face divergence B_n* A", e1, w);
            rec("0 harness = GIZMO run: B flux, face B states (per-face values)", DMAX(e2, e6), w);
            rec("0 harness = GIZMO run: DtB, Acc, DtU (running sums, scale = accumulated sum)", DMAX(DMAX(e3, e4), e5), w);
            double emax = DMAX(DMAX(DMAX(e1, e2), DMAX(e3, e4)), DMAX(e5, e6));
            emax = DMAX(DMAX(e1, e2), DMAX(DMAX(e3, e4), DMAX(e5, e6)));
            if(emax > 1e-10) {
                rec("0 count: faces with a reproduction difference > 1e-10 (listed as MISMATCH)", 1, w);
                printf("MISMATCH face=%ld step=%d eD=%.2e eFB=%.2e eDtB=%.2e eAcc=%.2e eDtU=%.2e eBstates=%.2e | P_M h/g %.17g %.17g | S_M h/g %.17g %.17g | BR0 h/g %.17g %.17g | DtU h/g %.17g %.17g\n",
                       id, r.ia[2], e1, e2, e3, e4, e5, e6, v[O_PM], vd[O_PM], v[O_SM], vd[O_SM], v[O_BR], vd[O_BR], v[O_DTU], vd[O_DTU]);
            }
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
            double e7 = DMAX(relerr_s(v + O_DTBI, vd + O_DTBI, vs + O_DTBI, 3), relerr_s(v + O_ADJI, vd + O_ADJI, vs + O_ADJI, 9));
            rec("0 harness = GIZMO run: psi adjoint DtB_i, gradB adjoint (running sums)", e7, w);
            emax = DMAX(emax, e7);
#endif
            if(emax > 1e-10) {exclude("row 0 failed (harness differs from GIZMO by > 1e-10): no derivative rows", w); continue;}
        }
        /* forward Jacobian of all outputs along the 24 directions */
        static double J[24][NOUT];
        for(int dir = 0; dir < 24; dir++)
        {
            LIn dL; PD dP[2]; CD dC[2]; memset(&dL, 0, sizeof(dL)); memset(dP, 0, sizeof(dP)); memset(dC, 0, sizeof(dC));
            seed(dir, &Cw[1], &dL, &dC[1]);
            double vv[NOUT], dv[NOUT]; memset(dv, 0, sizeof(dv));
            __enzyme_fwddiff<void>((void *)fvec, enzyme_dup, &r.L, &dL, enzyme_dup, Pw, dP, enzyme_dup, Cw, dC, enzyme_const, &aux, enzyme_dup, vv, dv);
            for(int q = 0; q < NOUT; q++) {J[dir][q] = dv[q];}
        }
        int dirk[3][8]; for(int k = 0; k < 3; k++) {dirk[k][0] = k; dirk[k][1] = 3 + k; for(int d = 0; d < 3; d++) {dirk[k][2 + d] = 6 + 3 * k + d; dirk[k][5 + d] = 15 + 3 * k + d;}}
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
        /* (b) hand-coded Jacobian of the face states (frozen wavespeed weights) vs Enzyme on the real reconstruction */
        double wL = v[O_WL], wR = v[O_WR];
        double hj[24]; for(int q = 0; q < 24; q++) {hj[q] = v[O_JAC + q];}
        rec("0 harness = GIZMO run: hand Jacobian psi_face_jac", relerr(hj, r.jac, 24), w);
        /* the face-state derivatives inside the face solve are not usable from Enzyme: Riemann_solver takes Riemann_vec by value and the
           shadow of that copy aliases the caller's (rotation and the HLLD B_n reset overwrite it). So differentiate the real reconstruction
           alone, with the face geometry of this face (hydro_core_meshless.h distance lines, copied) */
        double distance_from_i[3], distance_from_j[3];
        { struct kernel_hydra kernel = r.kh; double s_star_ij = 0, s_i, s_j, rinv = 1 / kernel.r; int k; /* s_star_ij = 0 as in hydro_core_meshless.h */
#include "ext_dist.h"
          (void)s_i; (void)s_j; }
        for(int k = 0; k < 3; k++)
        {
            double Bi = r.L.BPred[k], Bj = r.Cj0.Bfield()[k];
            double gi[3], gj[3]; for(int d = 0; d < 3; d++) {gi[d] = r.L.Gradients.B[k][d]; gj[d] = r.Cj0.Gradients.B[k][d];}
            double q1[2], q0[2]; recon_states(Bi, gi, Bj, gj, distance_from_i, distance_from_j, 1, q1); recon_states(Bi, gi, Bj, gj, distance_from_i, distance_from_j, 0, q0);
            double BLd = r.Rv.L.B[k], BRd = r.Rv.R.B[k];
            int second = (q1[0] == BLd && q1[1] == BRd), first = (q0[0] == BLd && q0[1] == BRd);
            double qlin_R = Bi + gi[0] * distance_from_i[0] + gi[1] * distance_from_i[1] + gi[2] * distance_from_i[2];
            double qlin_L = Bj + gj[0] * distance_from_j[0] + gj[1] * distance_from_j[1] + gj[2] * distance_from_j[2];
            int cat; const char *cn[] = {"unlimited", "limiter active", "first-order fallback", "equal cell values", "states not reproduced"};
            if(Bi == Bj) {cat = 3;} else if(second) {cat = (q1[0] == qlin_L && q1[1] == qlin_R) ? 0 : 1;} else if(first) {cat = 2;} else {cat = 4;}
            int mode = (cat == 2) ? 0 : 1;
            double je[8], jh[8];
            for(int m = 0; m < 8; m++)
            {
                double dQi = (m == 0), dQj = (m == 1), dgi[3] = {0, 0, 0}, dgj[3] = {0, 0, 0}, z3[3] = {0, 0, 0}, out2[2], dout2[2] = {0, 0};
                if(m >= 2 && m < 5) {dgi[m - 2] = 1;} if(m >= 5) {dgj[m - 5] = 1;}
                __enzyme_fwddiff<void>((void *)recon_states, Bi, dQi, enzyme_dup, gi, dgi, Bj, dQj, enzyme_dup, gj, dgj,
                                       enzyme_const, distance_from_i, enzyme_const, distance_from_j, enzyme_const, mode, enzyme_dup, out2, dout2);
                je[m] = 2 * wL * dout2[0] + 2 * wR * dout2[1]; jh[m] = hj[8 * k + m]; (void)z3;
            }
            {   /* the same Jacobian by Enzyme reverse mode */
                double qq[2] = {Bi, Bj}, bq[2] = {0, 0}, bgi[3] = {0, 0, 0}, bgj[3] = {0, 0, 0}, jr[8];
                __enzyme_autodiff<void>((void *)recon_comb, enzyme_dup, qq, bq, enzyme_dup, gi, bgi, enzyme_dup, gj, bgj, enzyme_const, distance_from_i,
                                        enzyme_const, distance_from_j, enzyme_const, mode, enzyme_const, 2 * wL, enzyme_const, 2 * wR);
                jr[0] = bq[0]; jr[1] = bq[1]; for(int d = 0; d < 3; d++) {jr[2 + d] = bgi[d]; jr[5 + d] = bgj[d];}
                if(cat != 4) {rec("(b) face-state Jacobian: Enzyme reverse = forward on reconstruct_face_states", relerr(jr, je, 8), w);}
            }
            cat_count[cat]++;
            if(cat != 0) {printf("BRANCH face=%ld step=%d k=%d %s\n", id, r.ia[2], k, cn[cat]);}
            char q[128]; snprintf(q, sizeof(q), "(b) face-state Jacobian hand vs Enzyme on reconstruct_face_states [%s]", cn[cat]);
            char w2[300]; snprintf(w2, sizeof(w2), "%s k=%d", w, k);
            if(cat != 4) {rec(q, relerr(jh, je, 8), w2);} else {rec("(b) EXCLUDED: face states not reproduced by the standalone reconstruction", 1, w2);}
            if(id == show) {printf("show k=%d cat=%s Bi=%.17g Bj=%.17g wL=%.17g wR=%.17g QL,QR=%.17g %.17g dump %.17g %.17g\n hand  ", k, cn[cat], Bi, Bj, wL, wR, q1[0], q1[1], BLd, BRd);
                            for(int m = 0; m < 8; m++) {printf(" %.17g", jh[m]);} printf("\n enzyme"); for(int m = 0; m < 8; m++) {printf(" %.17g", je[m]);} printf("\n");}
        }
#endif
        /* (a) adjoint: reverse mode on the centred face divergence (weights frozen by the product rule) vs the code's psi terms */
        double gD[24], gWL[24], gWR[24], gDfull[24];
        int outs[4] = {O_DC, O_WL, O_WR, O_D}; double *gs[4] = {gD, gWL, gWR, gDfull};
        for(int s = 0; s < 4; s++)
        {
#ifndef FWD_ONLY
            LIn bL; PD bP[2]; CD bC[2]; memset(&bL, 0, sizeof(bL)); memset(bP, 0, sizeof(bP)); memset(bC, 0, sizeof(bC));
            __enzyme_autodiff<void>((void *)fsel, enzyme_dup, &r.L, &bL, enzyme_dup, Pw, bP, enzyme_dup, Cw, bC, enzyme_const, &aux, enzyme_const, outs[s]);
            gather(&Cw[1], &bL, &bC[1], gs[s]);
#else
            for(int dir = 0; dir < 24; dir++) {gs[s][dir] = J[dir][outs[s]];} /* forward-mode columns (reverse mode did not compile) */
#endif
        }
#ifndef FWD_ONLY
        /* reverse vs forward consistency on D */
        { double gf[24]; for(int dir = 0; dir < 24; dir++) {gf[dir] = J[dir][O_D];} rec("Enzyme reverse = forward, d(face divergence)/d(B, gradB)", relerr(gDfull, gf, 24), w); }
#endif
        double x[24]; for(int k = 0; k < 3; k++) {x[k] = r.L.BPred[k]; x[3 + k] = r.Cj0.BPred[k] * r.Cj0.Density / r.Cj0.Mass;}
        for(int q = 0; q < 9; q++) {x[6 + q] = r.L.Gradients.B[q / 3][q % 3]; x[15 + q] = r.Cj0.Gradients.B[q / 3][q % 3];}
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
        {
            double fr[24], code[24], expct[24];
            for(int q = 0; q < 24; q++) {fr[q] = gD[q] + v[O_AN] * (v[O_BNL] * gWL[q] + v[O_BNR] * gWR[q]);}
            double psi_w = (r.L.PsiCons - (aux.jact ? H_Get_Gas_Psi_Conservative(Pw, Cw, 1) : 0.)) * All.cf_a2inv;
            for(int k = 0; k < 3; k++) {code[k] = v[O_DTBI + k]; code[3 + k] = aux.tbact ? v[O_DTBJ + k] : psi_w * fr[3 + k];}
            for(int q = 0; q < 9; q++) {code[6 + q] = v[O_ADJI + q]; code[15 + q] = aux.tbact ? v[O_ADJJ + q] : psi_w * fr[15 + q];}
            for(int q = 0; q < 24; q++) {expct[q] = psi_w * fr[q];}
            if(psi_w != 0) {rec("(a) psi/c_h adjoint: code DtB_PsiCons, PsiGradAdj vs psi * Enzyme dD_c/d(B,gradB)", relerr(code, expct, 24), w);}
            double s1 = 0, s1a = 0; for(int q = 0; q < 24; q++) {s1 += x[q] * code[q]; s1a += fabs(x[q] * code[q]);}
            double res = fabs(s1 - psi_w * v[O_DC]) / DMAX(s1a, fabs(psi_w * v[O_DC]));
            if(s1a > 0) {rec("(a) per-face centred identity sum x.adjoint = psi D_c (Euler)", res, w);}
            if(!aux.tbact && psi_w * v[O_DC] != 0)
            {   /* j inactive: the code adds no j-side adjoint (hydro_evaluate.h, psi_j_active); that share of psi D_c is not exchanged this step */
                double miss = 0; for(int k = 0; k < 3; k++) {miss += x[3 + k] * psi_w * fr[3 + k];} for(int q = 15; q < 24; q++) {miss += x[q] * psi_w * fr[q];}
                rec("(a) note: j inactive, j-side share of psi D_c not given to j: |share| / |psi D_c|", fabs(miss) / fabs(psi_w * v[O_DC]), w);
            }
            double ef = 0, es = 0; for(int q = 0; q < 24; q++) {ef = DMAX(ef, fabs(gD[q] - fr[q])); es = DMAX(es, fabs(fr[q]));}
            rec("(a) info: |unfrozen - frozen| / |frozen| (wavespeed-weight derivative, not used)", es > 0 ? ef / es : 0, w);
            if(id == show) {printf("show adjoint psi_w=%.17g D_c=%.17g\n", psi_w, v[O_DC]); for(int q = 0; q < 24; q++) {printf("  %2d x=%.17g code=%.17g enz=%.17g\n", q, x[q], code[q], expct[q]);}}
        }
#else
        {
            /* Dedner reference: B correction phi_mean A of this evaluation vs the pair adjoint (phi_i - phi_j) dD/dB_i of its own source */
            double phi_i = r.L.PhiPred * All.cf_a3inv, phi_j = H_Get_Gas_PhiField(Pw, Cw, 1) * All.cf_a3inv;
            double code[3], expct[3]; for(int k = 0; k < 3; k++) {code[k] = v[O_PHIM] * v[O_A + k]; expct[k] = (phi_i - phi_j) * gDfull[k];}
            if(phi_i != phi_j) {rec("(a) Dedner reference: phi_mean A vs (phi_i-phi_j) Enzyme dD/dB_i", relerr(code, expct, 3), w);}
        }
#endif
        /* (c) energy derivatives and bookings, per cell side, with the real energy function (estep_cell_mag_energy) */
        {
            PD Ps[2][2]; CD Cs[2][2];
            Ps[1][1] = Pw[1]; Cs[1][1] = Cw[1]; sync_cell(&Ps[1][1], &Cs[1][1]);
            cell_from_input(&r.L, &Pw[1], &Cw[1], &Ps[0][1], &Cs[0][1]); sync_cell(&Ps[0][1], &Cs[0][1]);
            for(int s = 0; s < 2; s++) {Ps[s][0] = Ps[s][1]; Cs[s][0] = Cs[s][1]; zero_acc(&Cs[s][1]); zero_acc(&Cs[s][0]);}
            double dEdB[2][3], dEdZ[2] = {0, 0}, psi_c[2] = {0, 0}, dZdD[2] = {0, 0};
            for(int s = 0; s < 2; s++)
            {
                PD bP[2]; CD bC[2]; memset(bP, 0, sizeof(bP)); memset(bC, 0, sizeof(bC));
                __enzyme_autodiff<void>((void *)Emag, enzyme_dup, Ps[s], bP, enzyme_dup, Cs[s], bC);
                Vec3<double> Bf = Cs[s][1].Bfield() * All.cf_a2inv; double eb[3], bb[3];
                for(int k = 0; k < 3; k++) {dEdB[s][k] = bC[1].B[k]; eb[k] = dEdB[s][k]; bb[k] = Bf[k];}
                rec("(c) dE_B/d(BV) (Enzyme, estep_cell_mag_energy) = Bfield() used in the u booking", relerr(eb, bb, 3), w);
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
                dEdZ[s] = bC[1].PsiZ;
                psi_c[s] = H_Get_Gas_Psi_Conservative(Ps[s], Cs[s], 1);
                /* dZ/dt per unit face divergence from the real final operations (divB = 1, nothing else) */
                PD Pt[2] = {Ps[s][1], Ps[s][1]}; CD Ct[2] = {Cs[s][1], Cs[s][1]}; zero_acc(&Ct[1]); Ct[1].divB = 1;
                fin_psi(Pt, Ct); dZdD[s] = Ct[1].DtPsiZ;
                double lhs = dEdZ[s] * dZdD[s], rhs = -psi_c[s] * All.cf_a2inv;
                rec("(c) dE_psi/dZ * dZ/dt per unit divB (Enzyme x code) = -psi used in the bookings", fabs(lhs - rhs) / DMAX(fabs(rhs), 1e-300), w);
                if(s == 0) {rec("(c) cell i rebuilt from hydro input: Get_Gas_Psi_Conservative = local.PsiCons", fabs(psi_c[0] - r.L.PsiCons) / DMAX(fabs(r.L.PsiCons), 1e-300), w);}
#endif
            }
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
            /* upwind dissipation of this face: Enzyme dE of both cells for the code's DtB_PsiDiss and PsiDivDiss (j mirrored), plus the heat booked */
            {
                double heat = 2 * (v[O_HEAT] - v[O_ADVHEAT]);
                double t[4];
                t[0] = 0; t[1] = 0; for(int k = 0; k < 3; k++) {t[0] += dEdB[0][k] * v[O_DTBDISS + k]; t[1] -= dEdB[1][k] * v[O_DTBDISS + k];}
                t[2] = dEdZ[0] * dZdD[0] * v[O_DIVDISS]; t[3] = dEdZ[1] * dZdD[1] * (-v[O_DIVDISS]);
                double dE = t[0] + t[1] + t[2] + t[3], sc = fabs(heat) + fabs(t[0]) + fabs(t[1]) + fabs(t[2]) + fabs(t[3]);
                if(sc > 0) {rec("(c) GLM upwind: Enzyme dE(B,psi) + heat booked = 0 per face", fabs(dE + heat) / sc, w);}
                rec("(c) GLM upwind heat >= 0 per face (value: -heat/scale if negative)", (heat < 0 && sc > 0) ? -heat / sc : 0, w);
                if(id == show) {printf("show upwind heat=%.17g dE terms %.17g %.17g %.17g %.17g\n", heat, t[0], t[1], t[2], t[3]);}
            }
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* MFV transport of Q = sqrt(m) Z: E_psi change of both cells (real final operations, Enzyme dE/dZ) plus the heat booked */
            if(v[O_MDOT] != 0)
            {
                double dz[2];
                for(int s = 0; s < 2; s++) {
                    PD Pt[2] = {Ps[s][1], Ps[s][1]}; CD Ct[2] = {Cs[s][1], Cs[s][1]}; zero_acc(&Ct[1]);
                    double sg = s ? -1 : 1; Ct[1].DtPsiQ = sg * v[O_DTPSIQ]; Ct[1].DtMass = sg * v[O_MDOT];
                    fin_psi(Pt, Ct); dz[s] = Ct[1].DtPsiZ;
                }
                double heat = 2 * v[O_ADVHEAT], dE = dEdZ[0] * dz[0] + dEdZ[1] * dz[1], sc = fabs(heat) + fabs(dEdZ[0] * dz[0]) + fabs(dEdZ[1] * dz[1]);
                if(sc > 0) {rec("(c) MFV psi transport: Enzyme dE_psi + heat booked = 0 per face", fabs(dE + heat) / sc, w);}
                rec("(c) MFV psi transport heat >= 0 per face (value: -heat/scale if negative)", (heat < 0 && sc > 0) ? -heat / sc : 0, w);
            }
#endif
#else
            /* Dedner: the centred phi correction phi_mean A on B_i is booked in u (ledger row 1); Enzyme dE_B along it */
            {
                double dE = 0; for(int k = 0; k < 3; k++) {dE += dEdB[0][k] * v[O_PHIM] * v[O_A + k];}
                double led = v[O_LEDG1]; double sc = DMAX(fabs(led), fabs(dE)), st = 0;
                for(int k = 0; k < 3; k++) {st += fabs(dEdB[0][k] * v[O_PHIM] * v[O_A + k]);}
                if(sc > 0) {rec("(c) Dedner phi_mean A: Enzyme dE_B = booked (ledger row 1)", fabs(dE - led) / sc, w);
                            rec("(c) Dedner phi_mean A: same, scale = sum of |terms| of the dot product", fabs(dE - led) / st, w);}
                if(id == show) {printf("show dedner-booking dE=%.17g ledger=%.17g sum|terms|=%.17g Bi(local)=%.17g %.17g %.17g dE/dBV=%.17g %.17g %.17g\n", dE, led, st,
                                       r.L.BPred[0], r.L.BPred[1], r.L.BPred[2], dEdB[0][0], dEdB[0][1], dEdB[0][2]);}
            }
#endif
            /* Powell terms of this face (cell i, divB = this face's B_n* A): Enzyme dE along the code's rates vs the ledger */
            {
                PD Pt[2] = {Ps[0][1], Ps[0][1]}; CD Ct[2] = {Cs[0][1], Cs[0][1]}; zero_acc(&Ct[1]); Ct[1].divB = v[O_D];
                fin_powell_convert(Pt, Ct);
                double du = Ct[1].DtInternalEnergy / Pt[1].Mass; /* the code's next line divides by the mass */
                PD dP[2]; CD dC[2]; memset(dP, 0, sizeof(dP)); memset(dC, 0, sizeof(dC));
                for(int k = 0; k < 3; k++) {dP[1].Vel[k] = Ct[1].HydroAccel[k]; dC[1].B[k] = Ct[1].DtB[k];}
                double ekb = __enzyme_fwddiff<double>((void *)Etot, enzyme_dup, Pt, dP, enzyme_dup, Ct, dC);
                dC[1].InternalEnergy = du;
                double etot = __enzyme_fwddiff<double>((void *)Etot, enzyme_dup, Pt, dP, enzyme_dup, Ct, dC);
#ifdef TEST_MHD_ENERGY_LEDGER_OUTPUT
                double led = Ct[1].DtE_MHDLedger[0];
                double sc = fabs(led) + fabs(ekb) + fabs(Pt[1].Mass * du);
                if(sc > 0) {
                    rec("(c) Powell: Enzyme d(E_kin+E_B+E_th) along code rates = ledger row 0", fabs(etot - led) / sc, w);
                    rec("(c) Powell: Enzyme d(E_kin+E_B) = ledger row 0 - m du (u booking)", fabs(ekb - (led - Pt[1].Mass * du)) / sc, w);
                }
                if(id == show) {printf("show powell D=%.17g ledger=%.17g dEkinB=%.17g m du=%.17g dEtot=%.17g\n", v[O_D], led, ekb, Pt[1].Mass * du, etot);}
#endif
                /* split by term: momentum (-B divB), induction (-v divB), thermal booking (total minus the two) */
                PD dPk[2]; CD dCk[2]; memset(dPk, 0, sizeof(dPk)); memset(dCk, 0, sizeof(dCk));
                for(int k = 0; k < 3; k++) {dPk[1].Vel[k] = Ct[1].HydroAccel[k];}
                double ek = __enzyme_fwddiff<double>((void *)Etot, enzyme_dup, Pt, dPk, enzyme_dup, Ct, dCk);
                memset(dPk, 0, sizeof(dPk)); for(int k = 0; k < 3; k++) {dCk[1].B[k] = Ct[1].DtB[k];}
                double eb = __enzyme_fwddiff<double>((void *)Etot, enzyme_dup, Pt, dPk, enzyme_dup, Ct, dCk);
                g_pw[0] += fabs(ek); g_pw[1] += fabs(eb); g_pw[2] += fabs(Pt[1].Mass * du); g_pw[3] += fabs(etot); g_pw[4] += ek; g_pw[5] += eb; g_pw[6] += Pt[1].Mass * du; g_pw[7] += etot;
                double tot = fabs(ek) + fabs(eb) + fabs(Pt[1].Mass * du);
                if(tot > 0) {rec("(P) Powell info: |kinetic work| / (|kin|+|mag|+|thermal|) per face", fabs(ek) / tot, w); rec("(P) Powell info: |net energy change| / (|kin|+|mag|+|thermal|) per face", fabs(etot) / tot, w);}
            }
            /* the Powell terms and the cleaning source read the same face divergence */
            rec("(P) Powell divB (out.divB contribution) = cleaning source face divergence (Fluxes.B_normal_corrected)", relerr(v + O_DIVB, v + O_D, 1), w);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            { double dpv = 0, vp = 0; for(int k = 0; k < 3; k++) {dpv = DMAX(dpv, fabs(r.Cj0.ParticleVel[k] - r.Cj0.VelPred[k])); vp = DMAX(vp, fabs(r.Cj0.VelPred[k]));}
              rec("(G) info MFV: |ParticleVel_j - VelPred_j| / |VelPred_j| (the pair uses VelPred for j, ParticleVel for i)", vp > 0 ? dpv / vp : dpv, w); }
#endif
            /* (G) Galilean invariance: boost V of all velocities, forward derivative at V = 0 (Enzyme), per axis k */
            {
                double R0[NR], T0[NT]; rates_i(&r.L, Pw, Cw, &aux, R0); powell_terms(&r.L, Pw, Cw, &aux, T0);
                double An = v[O_AN], D = v[O_D], Bi[3], b2 = 0; for(int k = 0; k < 3; k++) {Bi[k] = r.L.BPred[k] * All.cf_a2inv; b2 += Bi[k] * Bi[k];}
                double absB = sqrt(b2), Pi = r.L.Pressure * All.cf_a3inv, Vi = r.L.Mass / (r.L.Density * All.cf_a3inv), cs = DMAX(r.L.SoundSpeed, 1e-300);
                double Fsc = DMAX(An * (Pi + 0.5 * b2), gmax(v + O_FV, 3)); /* momentum-flux scale |A|(P + B^2/2) */
                double Bsc = fabs(D) + An * absB, DB = DMAX(fabs(D) * absB, 1e-300);
                for(int k = 0; k < 3; k++)
                {
                    LIn dL; PD dP[2]; CD dC[2]; double vv[NOUT], dv[NOUT], RR[NR], dR[NR], TT[NT], dT[NT];
                    boost_seed(k, &dL, dP, dC); memset(dv, 0, sizeof(dv));
                    __enzyme_fwddiff<void>((void *)fvec, enzyme_dup, &r.L, &dL, enzyme_dup, Pw, dP, enzyme_dup, Cw, dC, enzyme_const, &aux, enzyme_dup, vv, dv);
                    boost_seed(k, &dL, dP, dC); memset(dR, 0, sizeof(dR));
                    __enzyme_fwddiff<void>((void *)rates_i, enzyme_dup, &r.L, &dL, enzyme_dup, Pw, dP, enzyme_dup, Cw, dC, enzyme_const, &aux, enzyme_dup, RR, dR);
                    boost_seed(k, &dL, dP, dC); memset(dT, 0, sizeof(dT));
                    __enzyme_fwddiff<void>((void *)powell_terms, enzyme_dup, &r.L, &dL, enzyme_dup, Pw, dP, enzyme_dup, Cw, dC, enzyme_const, &aux, enzyme_dup, TT, dT);
                    char wk[300]; snprintf(wk, sizeof(wk), "%s k=%d", w, k);
                    /* face fluxes: exact Galilean laws */
                    double msc = DMAX(fabs(v[O_MDOT]), An * r.L.Density * All.cf_a3inv * cs);
                    rec("(G) mass flux: dF_m/dV = 0", fabs(dv[O_MDOT]) / msc, wk);
                    double rm[3]; for(int q = 0; q < 3; q++) {rm[q] = dv[O_FV + q] - ((q == k) ? v[O_MDOT] : 0);}
                    rec("(G) momentum flux: dF_p/dV_k - e_k F_m = 0", gmax(rm, 3) / Fsc, wk);
                    double rE = dv[O_FP] - v[O_FV + k];
                    rec("(G) energy flux: dF_E/dV_k - F_p,k = 0 (scale |A|(P+B^2/2))", fabs(rE) / Fsc, wk);
                    if(fabs(D) > 0) {rec("(G) info: energy-flux residual / (|D| |B_i|)", fabs(rE) / DB, wk);}
                    double rb[3]; for(int q = 0; q < 3; q++) {rb[q] = dv[O_FB + q] - ((q == k) ? D : 0);}
                    rec("(G) B flux: dF_B/dV_k - e_k D = 0 (the moving-face term -w B_n*)", gmax(rb, 3) / Bsc, wk);
                    rec("(G) face divergence D: dD/dV = 0", fabs(dv[O_D]) / Bsc, wk);
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
                    { int gi[] = {O_DTBI, O_DTBI + 1, O_DTBI + 2, O_DTBDISS, O_DTBDISS + 1, O_DTBDISS + 2, O_DIVDISS, O_HEAT, O_JUMP, O_WL, O_WR, O_DTPSIQ, O_ADVHEAT};
                      double dm = 0, sm = 0; for(int q : gi) {dm = DMAX(dm, fabs(dv[q])); sm = DMAX(sm, fabs(v[q]));} dm = DMAX(dm, gmax(dv + O_ADJI, 9)); sm = DMAX(sm, gmax(v + O_ADJI, 9));
                      rec("(G) psi/c_h rows (adjoint, gradB adjoint, upwind B/psi/heat, weights, MFV Q flux): d/dV = 0", sm > 0 ? dm / sm : dm, wk); }
#endif
#ifdef DIVBCLEANING_DEDNER_EVOLUTION
                    { int gi[] = {O_PHIM, O_PHIDB, O_PHICORR, O_PHICORR + 1, O_PHICORR + 2, O_LEDG1, O_DTPHI};
                      double dm = 0, sm = 0; for(int q : gi) {dm = DMAX(dm, fabs(dv[q])); sm = DMAX(sm, fabs(v[q]));}
                      rec("(G) Dedner rows (phi mean, B_n jump, PhiCorr, phi booking, MFV phi flux): d/dV = 0", sm > 0 ? dm / sm : dm, wk); }
#endif
                    /* Powell terms split */
                    rec("(G) Powell momentum -B divB: d/dV = 0", gmax(dT + T_MP, 3) / DB, wk);
                    double ri[3]; for(int q = 0; q < 3; q++) {ri[q] = dT[T_BP + q] + ((q == k) ? D : 0);}
                    rec("(G) Powell induction -v divB: d/dV_k = -e_k D exactly (cancels the B-flux term)", gmax(ri, 3) / Bsc, wk);
                    rec("(G) Powell energy -(v.B) divB: d/dV_k = -D B_k exactly", fabs(dT[T_EP] + D * Bi[k]) / DB, wk);
                    /* final cell rates (face share, real final operations) */
                    rec("(G) cell rate d(BV)/dt (flux + Powell + psi): d/dV = 0", gmax(dR + R_DB, 3) / Bsc, wk);
                    rec("(G) cell rate dv/dt: d/dV = 0", gmax(dR + R_DV, 3) / (Fsc / (r.L.Mass * cs)), wk);
                    rec("(G) cell rate m du/dt (incl. psi heat): d/dV = 0 (scale |A|(P+B^2/2))", fabs(dR[R_MDU]) / Fsc, wk);
#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH
                    rec("(G) cell rate dZ/dt: d/dV = 0", fabs(dR[R_DZ]) / (r.L.PsiCh * Bsc / sqrt(Vi) / cs), wk);
#endif
                    /* variants for the Powell energy */
                    rec("(V) code: d(m du)/dV / (|D||B_i|)", fabs(dR[R_MDU]) / DB, wk);
                    rec("(V) (a) Powell work into heat: d(m du + heat_a)/dV / (|D||B_i|)", fabs(dR[R_MDU] + dT[T_HA]) / DB, wk);
                    rec("(V) (a) info: d(heat_a)/dV_k / (2 D B_k)", fabs(Bi[k]) > 1e-3 * absB ? fabs(dT[T_HA] / (2 * D * Bi[k])) : 1, wk);
                    rec("(V) (b) Janhunen: m du_b - m du_code = -E_P + v.M_P (value / |D v B|)", fabs(T0[T_JAN]) / DMAX(fabs(D) * absB * DMAX(gmax(r.L.Vel.data, 3), 1e-300), 1e-300), wk);
                    rec("(V) (b) Janhunen: d(m du_b)/dV / (|D||B_i|)", fabs(dR[R_MDU] + dT[T_JAN]) / DB, wk);
                    rec("(V) (c) heat with v - w_face: d(heat_c)/dV / (|D||B_i|)", fabs(dT[T_HC]) / DB, wk);
                    rec("(V) (c) heat with v - w_face: d(m du + heat_c)/dV / (|D||B_i|)", fabs(dR[R_MDU] + dT[T_HC]) / DB, wk);
                }
                /* (c) booking: Enzyme d(E_kin + E_B + E_th) along the Powell rates + heat_c equals the row the code would book (v.M_P + B.dB_P + heat_c) */
                {
                    PD Pt[2]; CD Ct[2]; cell_from_input(&r.L, &Pw[1], &Cw[1], &Pt[1], &Ct[1]); sync_cell(&Pt[1], &Ct[1]); Pt[0] = Pt[1]; Ct[0] = Ct[1];
                    PD dP[2]; CD dC[2]; memset(dP, 0, sizeof(dP)); memset(dC, 0, sizeof(dC));
                    for(int k = 0; k < 3; k++) {dP[1].Vel[k] = T0[T_MP + k] / Pt[1].Mass; dC[1].B[k] = T0[T_BP + k];}
                    dC[1].InternalEnergy = T0[T_HC] / Pt[1].Mass;
                    double e = __enzyme_fwddiff<double>((void *)Etot, enzyme_dup, Pt, dP, enzyme_dup, Ct, dC);
                    double book = T0[T_W] + T0[T_HC], sc = fabs(T0[T_W]) + fabs(T0[T_HC]);
                    if(sc > 0) {rec("(V) (c) booking: Enzyme dE along Powell rates + heat_c = booked row (v.M_P + B.dB_P + heat_c)", fabs(e - book) / sc, w);}
                    (void)R0;
                }
            }
        }
    }
    printf("faces read %ld; face-state categories (per component): unlimited %ld, limited %ld, first-order %ld, equal %ld, not reproduced %ld\n",
           nf, cat_count[0], cat_count[1], cat_count[2], cat_count[3], cat_count[4]);
    printf("POWELL energy rates summed over the dumped faces (cell i side): sum|kin| %.6e sum|mag| %.6e sum|m du| %.6e sum|total| %.6e ; signed: kin %.6e mag %.6e thermal %.6e total %.6e\n",
           g_pw[0], g_pw[1], g_pw[2], g_pw[3], g_pw[4], g_pw[5], g_pw[6], g_pw[7]);
    int rc = cofactor_report("mhd_face");
    MPI_Finalize();
    return rc;
}
