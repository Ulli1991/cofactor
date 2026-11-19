"""mhd_face: the GIZMO hydro face (MFM/MFV, Dedner GLM or psi/c_h), its final per-cell operations and accessors.
Port of enzyme/mhd_audit/harness/gen.py (2026-11-18); the extracted text must be byte-identical to that script's output.
usage: python3 spec.py <build_dir> <out_dir> <manifest_dir> [--accept]"""
import sys

import extract  # cofactor lib/ on PYTHONPATH (the cofactor CLI sets it)

HD = "hydro/hydro_evaluate.h"
TL = "hydro/hydro_toplevel.cc"
FIN = "void hydro_final_operations_and_cleanup(void)"
SIG = "static double {name}(struct particle_data *Pw, struct gas_cell_data *Cw, int {arg})"


def spec(X, B):
    psi = "DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH" in (B / "Config.sh").read_text()
    # head of hydro_toplevel.cc: headers, statics fac_mu/fac_magnetic_pressure, Riemann structs, reimann.h, INPUT/OUTPUT structs
    X.block("ext_head.h", TL, "#include <mpi.h>", "static inline void particle2in_hydra(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration);")
    # hydro_force_evaluate: local declarations, per-i prelude, per-pair body (kernel.dp .. the dump hook TEST_FACE_DUMP_B)
    X.block("ext_decl.h", HD, "    int j, k, n, startnode, numngb, kernel_mode, listindex;", "    if(mode == 0)", start_after="int hydro_force_evaluate(")
    X.block("ext_prelude.h", HD, "    dt_hydrostep_i = local.dt_hydrostep_i;", "#ifdef RT_SOLVER_EXPLICIT", start_after="int hydro_force_evaluate(")
    X.block("ext_dist.h", "hydro/hydro_core_meshless.h", "s_i =  0.5 * kernel.r;",
            "for(k=0;k<3;k++) {distance_from_j[k] = distance_from_i[k] * s_j; distance_from_i[k] *= s_i;}", include_stop=True)
    X.block("ext_body.h", HD, "                kernel.dp = local.Pos - P[j].Pos;", "TEST_FACE_DUMP_B")
    # final operations, per cell (i): Powell terms; conversion of total to thermal energy; psi/c_h block
    X.block("ext_fin_powell.h", TL, "/* need to subtract out the source terms proportional", "double magnorm_closure = Get_DtB_FaceArea_Limiter",
            start_after=FIN)
    X.block("ext_fin_convert.h", TL, "/* we calculated the flux of conserved variables", "CellP[i].DtB *= magnorm_closure;", include_stop=True,
            start_after=FIN, expect_closed=1)
    if psi:
        X.block("ext_fin_psi.h", TL, "/* conjugate cleaning pair, added after the u booking",
                "{Vec3<double> Bl = CellP[i].Bfield() * All.cf_a2inv; CellP[i].DtE_MHDLedger[2]", start_after=FIN, expect_closed=1,
                cut_last="#ifdef TEST_MHD_ENERGY_LEDGER_OUTPUT", prefix="#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH\n", suffix="#endif\n")
    else:
        X.literal("ext_fin_psi.h", "\n")
    # accessors (core/predict.cc) and energies (core/energy_step_trace.cc), with P/CellP passed in
    PR, ES = "core/predict.cc", "core/energy_step_trace.cc"
    X.functions("ext_fun.h", [
        "#ifdef DIVBCLEANING_CONSERVATIVE_PSI_OVER_CH",
        (PR, r"^double Get_Gas_CleaningSpeed\(int i\)", SIG.replace("{name}", "H_Get_Gas_CleaningSpeed")),
        (PR, r"^double Get_Gas_Psi_Conservative\(int i\)", SIG.replace("{name}", "H_Get_Gas_Psi_Conservative")),
        (PR, r"^double INLINE_FUNC Get_Gas_PhiField\(int i\)", SIG.replace("{name}", "H_Get_Gas_PhiField")),
        "#else",
        (PR, r"^double INLINE_FUNC Get_Gas_PhiField\(int i_particle_id\)", SIG.replace("{name}", "H_Get_Gas_PhiField")),
        "#endif",
        (ES, r"^double estep_cell_mag_energy\(int i\)", SIG.replace("{name}", "H_estep_cell_mag_energy")),
        (ES, r"^double estep_cell_energy\(int i\)", SIG.replace("{name}", "H_estep_cell_energy")),
    ])


if __name__ == "__main__":
    extract.run(spec, sys.argv)
