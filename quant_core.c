/*quant_core.c*/
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <string.h>
#include "basic.h"
#include "index_core.h"
#include "quant_core.h"
#include "IO.h"
#include "cholmod_core.h"
#define BASISHS_DFLT (-1)

basis_Kref_t* build_basisKref (tref_cat_binary_t * btref, Kref_t *Kref){

    int nllong = btref->n % 32 == 0 ?  btref->n/32 :  (btref->n/32 + 1);
	int K = Kref->K;
	int Kref_size = Kref->size;
	int *ref = Kref->ref;

	int htsize = nextPrime( (int) ((double) Kref->size / LD_FCTR) );
	printf("basis hash table size = %d\n",htsize);

	//build basis_ht which map kmer to Kref	index
	int* basis_ht = malloc(sizeof (int) * htsize);
	for(int i = 0; i < htsize; i++) basis_ht[i] = BASISHS_DFLT ;
	
    for(int i = 0; i < Kref_size ;i++){
        int strnd  = ref[i] > 0 ? 1: -1;
        int abs_pos = abs(ref[i]) - 1;
        llong minkmer = pos2kmer(btref,nllong,strnd,abs_pos, K) ;

		for(int n = 0; n < htsize; n++){
            int hsv = HASH(minkmer,n,htsize);
			if( basis_ht[hsv] == BASISHS_DFLT){
				basis_ht[hsv] = i; // note here should be i, not ref[i] 
				break;
			}
		}	
    }		
	
	basis_Kref_t* basis_Kref = malloc(sizeof(basis_Kref_t)); 
	basis_Kref->htsize = htsize;
	basis_Kref->K = K; 
	basis_Kref->Kref_size = Kref_size;
	basis_Kref->basis_ht = basis_ht;
	basis_Kref->ref = ref;
	return (basis_Kref);
}

double * create_b (tref_cat_binary_t * btref, basis_Kref_t* basis_Kref, int fn, char *fq_f[]){
#define FQ_LEN 4096
int K = basis_Kref->K;
int comp_bittl = 64 - 2*K;
int crvsaddmove = 2*K - 2;
llong tupmask = _64MASK >> comp_bittl;
static llong tuple, crvstuple, unituple, base;
char seq[FQ_LEN];
char qual[FQ_LEN];
char ch;
FILE *fh;
int basenum; 

int bsize = basis_Kref->Kref_size;
double *b = calloc(bsize,sizeof(double));
int htsize = basis_Kref->htsize;
int *basis_ht = basis_Kref->basis_ht;
int *ref = basis_Kref->ref;
int nllong = btref->n % 32 == 0 ?  btref->n/32 :  (btref->n/32 + 1);
for(int i=0;i<fn;i++){
	fh = fopen(fq_f[i], "r");
	assert( (fh != NULL) && "File open failed in create_b()");
	while(fgets(seq,FQ_LEN,fh) && fgets(seq,FQ_LEN,fh) && fgets(qual,FQ_LEN,fh) && fgets(qual,FQ_LEN,fh)){
		base = 1;  
		for(int pos = 0; (ch = seq[pos]) != '\n'; pos++){
			basenum = Basemap[(int)ch];
			if(basenum != DEFAULT){
				tuple = ( ( tuple<< 2 ) | (llong)basenum ) & tupmask ;
				crvstuple = ( crvstuple >> 2 ) + (((llong)basenum^3LLU) << crvsaddmove);
				base++;	
			}else{ base = 1;continue;}			
			//this block is core: counting b
			if( base > K ){
				unituple = tuple < crvstuple ? tuple:crvstuple;	
				for(int n = 0; n < htsize; n++){
					int hsv = HASH(unituple,n,htsize);
					int index = basis_ht[hsv];
					if( index == BASISHS_DFLT ) break;
					else{
						int strnd = ref[index] > 0 ? 1: -1;
						int abs_pos = abs(ref[index]) - 1;
						llong minkmer = pos2kmer(btref,nllong,strnd,abs_pos, K) ;
						if(minkmer == unituple){
							b[index]++;
							break;
						}						
					}											
				}				
			}// core block end	
		}
	}// while 		
	fclose(fh);
}// for i	
return (b);

};

cholmod_dense * make_cholmod_dense (double * b, int nrow, cholmod_common *c){
	
	cholmod_dense * cholmod_b = malloc(sizeof(cholmod_dense));
	cholmod_b->nrow = nrow;
	cholmod_b->ncol = 1;
	cholmod_b->nzmax = nrow;
	cholmod_b->d = nrow;
	cholmod_b->x = b;
	cholmod_b->z = NULL;
	cholmod_b->xtype = CHOLMOD_REAL;
	return (cholmod_b);			
};























