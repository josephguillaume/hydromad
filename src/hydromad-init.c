#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* .C calls */
extern void ar1_tv(void *, void *, void *, void *, void *);
extern void filter_constloss(void *, void *, void *, void *, void *, void *);
extern void inverse_filter(void *, void *, void *, void *, void *, void *,
                           void *, void *, void *);
extern void inverse_filter_lambda(void *, void *, void *, void *, void *,
                                  void *, void *, void *, void *, void *);
extern void routing_gr4j(void *, void *, void *, void *, void *, void *, void *,
                         void *, void *);
extern void sma_awbm(void *, void *, void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *, void *, void *,
                     void *, void *);
extern void sma_bucket(void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *);
extern void sma_cmd(void *, void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *);
extern void sma_gr4j(void *, void *, void *, void *, void *, void *, void *,
                     void *);
extern void sma_sac(void *, void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *, void *, void *, void *);
extern void sma_sac_state(void *, void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *, void *,
                          void *, void *);
extern void sma_snow(void *, void *, void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *);
extern void sriv_system(void *, void *, void *, void *, void *, void *, void *,
                        void *, void *, void *);
extern void swimp_core(void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _hydromad_simhyd_sim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                 SEXP, SEXP, SEXP);

extern SEXP _hydromad_hbv_sim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                              SEXP, SEXP, SEXP, SEXP);

extern SEXP _hydromad_hbvrouting_sim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                     SEXP, SEXP);
extern SEXP _hydromad_hbv_pet(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"ar1_tv", (DL_FUNC)&ar1_tv, 5},
    {"filter_constloss", (DL_FUNC)&filter_constloss, 6},
    {"inverse_filter", (DL_FUNC)&inverse_filter, 9},
    {"inverse_filter_lambda", (DL_FUNC)&inverse_filter_lambda, 10},
    {"routing_gr4j", (DL_FUNC)&routing_gr4j, 9},
    {"sma_awbm", (DL_FUNC)&sma_awbm, 16},
    {"sma_bucket", (DL_FUNC)&sma_bucket, 13},
    {"sma_cmd", (DL_FUNC)&sma_cmd, 11},
    {"sma_gr4j", (DL_FUNC)&sma_gr4j, 8},
    {"sma_sac", (DL_FUNC)&sma_sac, 14},
    {"sma_sac_state", (DL_FUNC)&sma_sac_state, 32},
    {"sma_snow", (DL_FUNC)&sma_snow, 12},
    {"sriv_system", (DL_FUNC)&sriv_system, 10},
    {"swimp_core", (DL_FUNC)&swimp_core, 19},
    {NULL, NULL, 0}};

static const R_CallMethodDef CallEntries[] = {
    {"_hydromad_simhyd_sim", (DL_FUNC)&_hydromad_simhyd_sim, 11},
    {"_hydromad_hbv_sim", (DL_FUNC)&_hydromad_hbv_sim, 12},
    {"_hydromad_hbvrouting_sim", (DL_FUNC)&_hydromad_hbvrouting_sim, 9},
    {"_hydromad_hbv_pet", (DL_FUNC)&_hydromad_hbv_pet, 5},
    {NULL, NULL, 0}};

void R_init_hydromad(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
