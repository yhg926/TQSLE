#include "index_core.h"
typedef struct basisKref {
    int htsize;
    int K;
    int Kref_size;
    int *basis_ht;
    int *ref;
} basis_Kref_t ;

basis_Kref_t* build_basisKref (tref_cat_binary_t *, Kref_t *);
double* create_b (tref_cat_binary_t*, basis_Kref_t*, int, int, char **);
double * create_b2 (pos2bidx_t *, tref_cat_binary_t *, basis_Kref_t*, int, int, char **);
cholmod_dense * make_cholmod_dense (double*, int, cholmod_common *);
