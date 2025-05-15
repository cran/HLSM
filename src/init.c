#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void sampleFixedIntervention(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sampleRandomIntervention(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sampleFixedEFFixedIntercept(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sampleFixedIntercept(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"sampleFixedIntervention",  (DL_FUNC) &sampleFixedIntervention,  28},
  {"sampleRandomIntervention", (DL_FUNC) &sampleRandomIntervention, 30},
    {"sampleFixedEFFixedIntercept", (DL_FUNC) &sampleFixedEFFixedIntercept, 25},
    {"sampleFixedIntercept", (DL_FUNC) &sampleFixedIntercept, 25},
  {NULL, NULL, 0}
};

void R_init_HLSM(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
