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

double * create_b2 (pos2bidx_t * Pos2bidx, tref_cat_binary_t * btref, basis_Kref_t* basis_Kref, int fn, int is_fa, char *fq_f[]){
#define THREAD_MAX 65536
#define FQ_LEN 4096
#define BFQ_LEN (FQ_LEN/32) //for llong type
	char (*fq_buff)[FQ_LEN] = malloc( THREAD_MAX * FQ_LEN );
	char tmp[FQ_LEN];

	llong (*bfq_buff)[BFQ_LEN] = malloc(THREAD_MAX * BFQ_LEN *sizeof(llong));
//int len = Pos2bidx->len;
	int *pos2bidx = Pos2bidx->pos2bidx;
	int K = basis_Kref->K;
	int lshft = (32-K)*2;
	//int comp_bittl = 64 - 2*K;
	//int crvsaddmove = 2*K - 2;
	//llong tupmask = _64MASK >> comp_bittl;
	FILE *fh;
	int bsize = basis_Kref->Kref_size;
	double *b = calloc(bsize,sizeof(double));
	int htsize = basis_Kref->htsize;
	int *basis_ht = basis_Kref->basis_ht;
	int *ref = basis_Kref->ref;
	int nllong = btref->n % 32 == 0 ?  btref->n/32 :  (btref->n/32 + 1);
	int l;

	llong *cat_tref = btref->cat_tref;
	llong *rc_cat_tref = btref->rc_cat_tref; 

	if(is_fa)
  	printf ("[-f %d] The RNA-seq data is assumed to be fasta format\n",is_fa);


	for(int fi=0;fi<fn;fi++){
  	fh = fopen(fq_f[fi], "r");
  	assert( (fh != NULL) && "File open failed in create_b()");

  	while (!feof(fh)) {
    	if(is_fa)
      	for (l = 0; l < THREAD_MAX && fgets(tmp,FQ_LEN,fh) && fgets(fq_buff[l],FQ_LEN,fh) ; l++) ;
    	else
      	for (l = 0; l < THREAD_MAX && fgets(tmp,FQ_LEN,fh) && fgets(fq_buff[l],FQ_LEN,fh)
					 && fgets(tmp,FQ_LEN,fh) && fgets(tmp,FQ_LEN,fh); l++) ;

#pragma omp parallel for schedule(guided)
    	for (int t = 0 ; t < l ; t++ ){
				int rlen = strlen(fq_buff[t]) - 1; // \n is omitted
				int rlenllong = rlen % 32 == 0 ? rlen/32 : (rlen/32 + 1)	;		

				//create binary reads/
				llong temp = 0LLU;
				int j;
				for (j = 0; j < rlen; j++ ){
					temp <<=2;
					int basenum = Basemap[(int)fq_buff[t][j]] ;
				//caution: set non ACGT character to 0(A) 
					if(basenum != DEFAULT) 
						temp += (llong)basenum;

					if( j % 32 == 31 ) 
						bfq_buff[t][(int)(j / 32)]	= temp;				
				}      		
				if( j % 32 < 31 ) 
					bfq_buff[t][(int)(j / 32)]  = ( temp << 2*(32 - (j % 32) )) ;  


				int pos = 0;
				while( pos < rlen - K ){
					int match_len = 0; //read to tfref match length
					int ind = pos/32 ;
					int rmd = (pos % 32)*2;	
					llong tuple = (bfq_buff[t][ind]  << rmd ) >>lshft  ;
					if(rmd > lshft)
						tuple += bfq_buff[t][ind+1] >> (64 - rmd + lshft) ;

					llong crvstuple = crvs64bits(tuple) >> lshft;

					int read_strnd ; //1 for +, -1 for -
					llong unituple;

					if (tuple > crvstuple) {
         		read_strnd = -1 ;
          	unituple = crvstuple;
        	}
					else{
						read_strnd = 1 ;
						unituple = tuple;
					}
				
					for(int n = 0; n < htsize; n++){
						//int hsv = HASH(unituple,n,htsize);
						int index = basis_ht[HASH(unituple,n,htsize)];// int index = basis_ht[hsv];
						if( index == BASISHS_DFLT ) break;
						else { 
							int strnd = ref[index] > 0 ? 1: -1;
							int abs_pos = abs(ref[index]) - 1;
							int tind = abs_pos / 32 ;
							int trmd = (abs_pos % 32) * 2; //btref index, reminder
							llong minkmer;
				
							if(strnd == 1) {//+ strand
								minkmer = cat_tref[tind] << trmd >> lshft ;
								if( trmd > lshft )
									minkmer += cat_tref[tind+1] >> (64 - trmd + lshft);		
							}
							else {          //- strand
								minkmer = ( rc_cat_tref[nllong - 1 - tind] >>trmd) & (_64MASK >> lshft) ;
								if( trmd > lshft )
									minkmer += (rc_cat_tref[nllong - 2 - tind] << (64-trmd+lshft) ) >> lshft;				
							}
																								
							if(minkmer == unituple){
								//updata ind and rmd so as to skip the matched kmer						
								ind += (int)((rmd + 2*K)/64);
								rmd = (rmd + 2*K) % 64;																						
							
								llong *uni_tref;							
								if(read_strnd == strnd) { //same strand
									tind += (int)((trmd + 2*K)/64);
              		trmd = (trmd + 2*K) % 64;	
									uni_tref = cat_tref;						
								}
								else{ // use rc_cat_tref strand
									tind = nllong - 1 - tind;
									trmd = 64 - trmd;
                	uni_tref = rc_cat_tref;								
								}
							
								int cmp_nllong = (rlenllong - ind) < (nllong - tind) ? rlenllong - ind : nllong - tind;	
								//find match length / 							
								for (int i = 0; i< cmp_nllong ;i++){
									llong r_xo_f = ( (bfq_buff[t][ind+i] << rmd ) | ( bfq_buff[t][ind+i+1] >> (64-rmd) ) )    								
									       ^ ( (uni_tref[tind + i] << trmd) | (uni_tref[tind + i + 1] >> (64-trmd)  ) );

									if(r_xo_f) {
										int tmp_match_len = 16;									
										for(int bsft = 8; bsft > 0 ; bsft >>= 1  ){

											if(  r_xo_f >> (64 - tmp_match_len*2 ) )
												tmp_match_len -= bsft;					
											else
												tmp_match_len += bsft;			 				
 	
										}

										if(  r_xo_f >> (64 - tmp_match_len*2 ) ) 
											tmp_match_len -=1;
										
									 	match_len += tmp_match_len ;
																						
										break;
									} // r_xo_f != 0
									else
										 match_len += 32;								 
								} // find match length					
								// count b	
							//	printf("abs_pos=%d\tread_strnd=%d\tstrnd=%d\tmatch_len=%d\tpos=%d\n",abs_pos,read_strnd,strnd,match_len,pos);									
								if( read_strnd ==  strnd){
									for (int p = 0; p < match_len + 1 && pos + p < rlen -K && abs_pos < btref->n ; p++ ){
#pragma omp atomic
            				b[pos2bidx[abs_pos++]]++;
          				}
								}
								else{
									for (int p = 0; p < match_len + 1 && pos + p < rlen -K && abs_pos >0 ; p++ ){
#pragma omp atomic
                    b[pos2bidx[abs_pos--]]++;
                  }
								}
								pos += match_len;				
								break;		
							} // hash found element							
						} // hash collsion						
					} // for hashing  

					pos++ ;		
				} //while read pos
			}// for threads 

		}// while reads buffer

		fclose(fh);		
	}// for files loop

	free(fq_buff);
	free(bfq_buff);
	return (b);
};


/*

inline void b_kmer_count ( char **fq_buff, double *b) {
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
                              b[index]++;
                              break;
                          }
                      }
                  }
              }// core block end
          }//for
}


*/










































