/*quant_core.c*/
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
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

double * create_b (tref_cat_binary_t * btref, basis_Kref_t* basis_Kref, int fn, int is_fa, char *fq_f[]){
#define THREAD_MAX 65536
#define FQ_LEN 4096
char (*fq_buff)[FQ_LEN] = malloc( THREAD_MAX * FQ_LEN );
char tmp[FQ_LEN];

int K = basis_Kref->K;
int comp_bittl = 64 - 2*K;
int crvsaddmove = 2*K - 2;
llong tupmask = _64MASK >> comp_bittl;
FILE *fh;

int bsize = basis_Kref->Kref_size;
double *b = calloc(bsize,sizeof(double));
int htsize = basis_Kref->htsize;
int *basis_ht = basis_Kref->basis_ht;
int *ref = basis_Kref->ref;
int nllong = btref->n % 32 == 0 ?  btref->n/32 :  (btref->n/32 + 1);
int l;

if(is_fa)
	printf ("[-f %d] The RNA-seq data is assumed to be fasta format\n",is_fa);	


for(int i=0;i<fn;i++){
	fh = fopen(fq_f[i], "r");
	assert( (fh != NULL) && "File open failed in create_b()");

	while (!feof(fh)) {
		if(is_fa)
			for (l = 0; l < THREAD_MAX && fgets(tmp,FQ_LEN,fh) && fgets(fq_buff[l],FQ_LEN,fh) ; l++) ;
		else
			for (l = 0; l < THREAD_MAX && fgets(tmp,FQ_LEN,fh) && fgets(fq_buff[l],FQ_LEN,fh) && fgets(tmp,FQ_LEN,fh) && fgets(tmp,FQ_LEN,fh); l++) ; 
		
#pragma omp parallel for schedule(guided) 
		for (int t = 0 ; t < l ; t++ ){
			int base = 1; char ch; 
			llong tuple = 0LLU;
			llong crvstuple = 0LLU;
			llong unituple = 0LLU;
        	for(int pos = 0; (ch = fq_buff[t][pos]) != '\n'; pos++){
            	int basenum = Basemap[(int)ch];
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
#pragma omp atomic //#pragma omp critical is too slow
                            	b[index]++;
                            	break;
                        	}
                    	}
                	}
            	}// core block end
        	}//for								
		}//parallel for		
	}//while
	fclose(fh);
}// for i	

free(fq_buff);
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























