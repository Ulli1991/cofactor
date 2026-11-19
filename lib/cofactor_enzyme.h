// cofactor: Enzyme entry points (forward and reverse mode) and activity markers; no GIZMO dependence.
#ifndef COFACTOR_ENZYME_H
#define COFACTOR_ENZYME_H
extern int enzyme_dup, enzyme_const, enzyme_dupnoneed, enzyme_out;
template <typename R, typename... A> R __enzyme_fwddiff(void *, A...);
template <typename R, typename... A> R __enzyme_autodiff(void *, A...);
#endif
