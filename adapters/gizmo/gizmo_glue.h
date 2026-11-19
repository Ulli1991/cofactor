// cofactor adapter for GIZMO: stubs and inactive registrations for a harness that compiles real GIZMO code (needs allvars.h, MPI).
// Include after <cstdio> and the GIZMO declarations (allvars.h), before any extracted code.
#ifndef COFACTOR_GIZMO_GLUE_H
#define COFACTOR_GIZMO_GLUE_H
#include <cstdio>
#include <cstdlib>

/* GIZMO's warning/abort paths: counted, never executed; the record that reaches one is excluded by the harness */
static int g_warn_hits = 0, g_abort_hits = 0;
#undef fflush
#define fflush(x) (g_warn_hits++, 0) /* GIZMO calls fflush only after printing an error/warning */
#undef PRINT_WARNING
#define PRINT_WARNING(...) {g_warn_hits++;}
#undef endrun
#define endrun(x) {g_abort_hits++;}

#include "cofactor_enzyme.h"
/* globals and library calls without derivatives */
__attribute__((used)) void *__enzyme_inactive_global_All = (void *)&All;
__attribute__((used)) void *__enzyme_inactivefn_1 = (void *)&fflush;
__attribute__((used)) void *__enzyme_inactivefn_2 = (void *)&printf;
__attribute__((used)) void *__enzyme_inactivefn_3 = (void *)&snprintf;
__attribute__((used)) void *__enzyme_inactivefn_4 = (void *)&MPI_Abort;
__attribute__((used)) void *__enzyme_inactivefn_5 = (void *)&exit;

typedef struct particle_data PD;
typedef struct gas_cell_data CD;
#endif
