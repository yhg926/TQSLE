#ifndef INDEX_H
#define INDEX_H

#include "basic.h"
#include "cholmod.h"
#include "iqsle.h"

typedef struct trscrpt_ref {
    char * name;
    char * seq;
    int len ;
} trscrpt_ref_t ;

typedef struct trscrpt_ref_arr {
	int n; //# of ref transcript
	trscrpt_ref_t * trscrpt_ref_idx;
} trscrpt_ref_arr_t ; 
/*moved to basic.h 
typedef struct tref_cat_binary {
	int n; //length of catted ref ( # of nt )
	llong *cat_tref;
	llong *rc_cat_tref; //reverse complete
	int noN; //numer of invalid bases;
	int * Npos; //invalid base positions
	int noT; // numer of tref;
	int *tref_nxt_start; // i+1 trnscrpt start pos, can fast caculate i_start_pos

} tref_cat_binary_t ;
*/

//#define FST_POS_DFLT 0
#define KMER_HSH_RNK_DFLT (-1)
#define POS_INVLD (-2)

typedef struct btref_kmer_index {
	int tid; // which trnscrpt this position belong to;
	int kmer_hash_rank;//start from 0, -1 if kmer if uniq in tref; -2 if this pos is invalid e.g. contain N..., default = -1
} btref_kmer_index_t;

typedef struct ext_btref_kmer_index{
	int K;
	int htsize;
	int uniq_kmer_n;
	int npos; //# of valid pos;
	btref_kmer_index_t * btref_kmer_index;
	int *kmer_pos_ht;
	int **kmer_pos_mtx;// include npos and *pos; //all pos that shared this kmer, element take minus if minor kmer from minus
} ext_btref_kmer_index_t; // btref_kmer_index with a kmer hash table and a non-redundant kmer-position shared array

typedef struct Twght {
	int tid;
	int wght; //weight sum
} Twght_t;

typedef struct kmer_Twght_mtx{
    int nzmax;
    int nrow;
    int ncol;
    int *p; //row start pointers, size = nrow + 1 
    int *i; //col indices(tids), size = nzmax
    double *x; //weights, size = nzmax
} kmer_Twght_mtx_t;

typedef struct Kref {
	int K;
	int size;
	int *ref;//definition of row id, with size of nrow
} Kref_t;

typedef struct Kref_mtx{
	kmer_Twght_mtx_t * kmer_Twght_mtx;
	Kref_t *Kref; 
} Kref_mtx_t;


trscrpt_ref_arr_t *Read_fas_trscrpt (char *, int );

void trscrpt_ref_arr_free ( trscrpt_ref_arr_t * );

tref_cat_binary_t * cat_tref_arr(trscrpt_ref_arr_t * ) ;
void crvs_btref(tref_cat_binary_t * btref);

//set btref_kmer_index.tid, set .kmer_hash_rank to -1 if invalid set to -2
btref_kmer_index_t* build_btref_kmer_index ( tref_cat_binary_t*, int);
ext_btref_kmer_index_t * build_ext_btref_kmer_index (btref_kmer_index_t*, tref_cat_binary_t*, int) ;
Kref_mtx_t* build_Kref_mtx (ext_btref_kmer_index_t*, tref_cat_binary_t*, indexProperty_t*);
cholmod_sparse* kmer_Twght_mtx2cholmod_sparse(kmer_Twght_mtx_t*);

void free_btref(tref_cat_binary_t * );


#endif





