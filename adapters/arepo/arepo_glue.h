// cofactor adapter for AREPO: abort/warning/print paths of copied AREPO code are counted, never executed.
// Include before the extracted code. The record (call) that reaches one of them is excluded or reported by the harness.
#ifndef COFACTOR_AREPO_GLUE_H
#define COFACTOR_AREPO_GLUE_H
static int g_terminate_hits = 0, g_warn_hits = 0;
#undef Terminate
#define Terminate(...) (g_terminate_hits++)
#undef warn
#define warn(...) (g_warn_hits++)
namespace arepo
{
template <typename... A> static inline void print(A...) {}
}
#endif
