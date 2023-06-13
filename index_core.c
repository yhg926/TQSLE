#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "basic.h"
#include "cholmod.h"
#include "IO.h"

#define MAX_TRANSCRIPTS_LEN 500000
#define MAX_Ref_NUM 1000000

trscrpt_ref_arr_t *Read_fas_trscrpt (char *path, int min_len ){
	int c,refnum,refindex = 0, nskip = 0;
	FILE *fs;
	if ( NULL == (fs = fopen(path,"r")) )  perror("Error in opening file");

	trscrpt_ref_arr_t *trscrpt_ref_arr = malloc(sizeof(trscrpt_ref_arr_t)) ;

	char *line = malloc(MAX_TRANSCRIPTS_LEN);
	
	trscrpt_ref_arr->trscrpt_ref_idx = malloc( sizeof(trscrpt_ref_t) * MAX_Ref_NUM );
	

	while ( ( (c = fgetc(fs)) != '>') && (c != EOF)  );

	for( refnum = 0; c != EOF ; refindex++ ){
		if (refnum > MAX_Ref_NUM) { printf("# references exceed MAX_Ref_NUM %d: aborted!\n", MAX_Ref_NUM); exit(1); }

		if( (c=='>') && (c != EOF) ) {
			int pos_name;
			for ( pos_name = 0 ; ( (c = fgetc(fs)) != '\n') && (c != EOF) ; pos_name++ ){

				if(pos_name >= MAX_TRANSCRIPTS_LEN){
					line[ MAX_TRANSCRIPTS_LEN - 1 ] = '\0';	
                    printf("%dth transcript: %s exceed MAX_TRANSCRIPTS_LEN %d: aborted\n", refindex + 1, line, MAX_TRANSCRIPTS_LEN);
                    exit(1);
                }				
				line[pos_name] = c ;		
			};
			if (c == '\n') {
				(trscrpt_ref_arr->trscrpt_ref_idx)[refnum].name = realloc ( (trscrpt_ref_arr->trscrpt_ref_idx)[refnum].name, pos_name + 1 );
				char *temp = strchr(line,' '); 
				if (temp != NULL) {*temp = '\0';}
				else{
					temp = strchr(line,'|');
					if (temp != NULL) {*temp = '\0';}
				}
				strncpy((trscrpt_ref_arr->trscrpt_ref_idx)[refnum].name, line, pos_name);		
			} 
			for ( pos_name = 0 ; ( (c = fgetc(fs)) != '>') && (c != EOF) ; ){
				if(pos_name > MAX_TRANSCRIPTS_LEN){
					 printf("%dth transcript %s exceed MAX_TRANSCRIPTS_LEN %d: aborted\n", refindex + 1,
						(trscrpt_ref_arr->trscrpt_ref_idx)[refnum].name, MAX_TRANSCRIPTS_LEN);
					 exit(1);
				}					
				if ( (c!='\n')){						
					line[pos_name] = c ;
					pos_name++;
				}
			}
			if (pos_name < min_len) {
				nskip++;
				//printf("skipped %dth transcript: %s due to length %d is shorter than MIN_LEN %d\n", refindex + 1, (trscrpt_ref_arr->trscrpt_ref_idx)[refnum].name, pos_name, min_len);
				continue;
			}

			(trscrpt_ref_arr->trscrpt_ref_idx)[refnum].seq =  malloc(pos_name + 1);
			strncpy((trscrpt_ref_arr->trscrpt_ref_idx)[refnum].seq, line, pos_name);
			((trscrpt_ref_arr->trscrpt_ref_idx)[refnum].len) = pos_name;
			refnum++;	
		};
	}

	trscrpt_ref_arr->n = refnum;
	free(line);

	printf("total %d references found, %d is effective, %d is skipped due to shorter than fragment length (%d)\n",refindex,refnum,nskip, min_len);

	return(trscrpt_ref_arr); 			
}

void trscrpt_ref_arr_free ( trscrpt_ref_arr_t * tref ){

	for(int i = 0; i< tref->n; i++){
		free((tref->trscrpt_ref_idx)[i].name);
		free((tref->trscrpt_ref_idx)[i].seq);
	}	

	free(tref->trscrpt_ref_idx);
	free(tref);
}

void crvs_btref(tref_cat_binary_t * btref) { //creat cmplt rvs catted bref if cmplt rvs strand is NULL
	assert( (btref->n != 0) && (btref->cat_tref != NULL) ); 

	if(btref->rc_cat_tref == NULL){
		int arr_size; //size of catted_tref->rc_cat_tref in # of llong element
		if( btref->n % 32 == 0) 
			arr_size = btref->n / 32;
		else
			arr_size = (int) (btref->n / 32) + 1;
	
		btref->rc_cat_tref =  (llong *)malloc( arr_size*sizeof(llong) )  ;
		for(int i = 0; i < arr_size; i ++ )
        	(btref->rc_cat_tref)[arr_size - 1 - i] = crvs64bits((btref->cat_tref)[i]);
	}
		
}

tref_cat_binary_t * cat_tref_arr(trscrpt_ref_arr_t * tref_arr){

	tref_cat_binary_t * btref = malloc(sizeof(tref_cat_binary_t));
	/*initialize*/
	btref->n = 0;
	btref->cat_tref = NULL;
	btref->rc_cat_tref = NULL;
	btref->noN = 0;
	btref->Npos = NULL;	

	btref->noT = tref_arr->n; //get # of tref	
	btref->tref_nxt_start = malloc(sizeof(int) * (btref->noT) ); //index next tref start pos

	llong temp = 0LLU; // current 32-mer binary
	int pos = 0;		// current 32-mer start pos
	int base = DEFAULT; // current base	

	for(int i = 0; i < btref->noT; i++){
		int len = (tref_arr->trscrpt_ref_idx)[i].len;
		int tref_i_start_pos = btref->n ; // tref[i]'s nt start pos in the catted nt sequence

		btref->n += len; 
		(btref->tref_nxt_start)[i] = btref->n ; // record tref[i+1]'s nt start pos in the catted nt sequence

	 	btref->cat_tref = (llong *)realloc(btref->cat_tref,  ((int)(btref->n /32) + 1)*sizeof(llong) ); //8bits -> 2 bits
				

		char *seq = (tref_arr->trscrpt_ref_idx)[i].seq;

		for (int j = 0; j < len; j++ ){
			pos = tref_i_start_pos + j;
			temp <<=2;		
			base = Basemap[(int)seq[j]] ;
			if( base != DEFAULT )	temp += (llong)base; 
			else {
				btref->Npos = realloc(btref->Npos, (int)( btref->noN / 1024 + 1) * 1024 *sizeof(int)  );
				(btref->Npos)[btref->noN] = pos;
				btref->noN++;
			}	
		
			if( pos % 32 == 31 ) // make sure 32*2 is bit size of llong, the type of catted_tref->cat_tref element; 
				(btref->cat_tref) [(int)(pos / 32)] = temp;
		}
	}
	int rmd = pos % 32; // take reminder;
	//20230613:fixed bug:64-2*rmd -> 64-2*(rmd+1)
	if( rmd < 31 ) (btref->cat_tref)[(int)(pos / 32) ] = ( temp << (64-2*(rmd+1))) ;
	/*caculate rev_complete binary*/	
	crvs_btref(btref);		

	return(btref);
}

void free_btref(tref_cat_binary_t * btref){
	free(btref->cat_tref);
	free(btref->rc_cat_tref);
	free(btref->Npos);
	free(btref->tref_nxt_start);
	free(btref);
}

//malloc mem for btref_kmer_index and initialize set btref_kmer_index.tid
btref_kmer_index_t* build_btref_kmer_index ( tref_cat_binary_t *btref, int K){

	btref_kmer_index_t* btref_kmer_index = malloc( btref->n * sizeof(btref_kmer_index_t) );	
/*initialize btref_kmer_index*/
/* set btref_kmer_index.tid, set .kmer_hash_rank to -1 if invalid set to -2*/
	for(int i = 0 ; i< btref->noT; i++){
		int ith_tref_start = i == 0? 0:(btref->tref_nxt_start)[i-1];
		for(int j = ith_tref_start; j< (btref->tref_nxt_start)[i] ; j++ ){			
			btref_kmer_index[j].tid = i;
			btref_kmer_index[j].kmer_hash_rank = KMER_HSH_RNK_DFLT ;//-1
		}
	}		
//disable postions K-1 ahead '[^ACGT]' base
	for(int n = 0 ; n< btref->noN;n++){
		int Npos = (btref->Npos)[n] ;
		int tid = btref_kmer_index[Npos].tid ;
		int tstart;
		if(tid == 0) tstart = 0;
		else tstart =  (btref->tref_nxt_start)[tid - 1];

		int Nstart;
		if ( Npos - K + 1 < tstart) Nstart = tstart;	
		else Nstart = Npos - K + 1 ;

		for (int i = Nstart ; i <= Npos ; i++)
			btref_kmer_index[i].kmer_hash_rank = POS_INVLD;		
		
	}
//disable postions K-1 ahead tref i+1 start
	for(int n = 0 ; n< btref->noT;n++){
        int Npos = (btref->tref_nxt_start)[n] ;
        int tstart;
        if(n == 0) tstart = 0;
        else tstart =  (btref->tref_nxt_start)[n - 1];

        int Nstart;
        if ( Npos - tstart < K ) Nstart = tstart;
        else Nstart = Npos - K + 1 ;

        for (int i = Nstart ; i <  Npos ; i++)
            btref_kmer_index[i].kmer_hash_rank = POS_INVLD;

    }
	
	return (btref_kmer_index);
}
/* has moved to basic.h (needed by quant_core.c)
static inline llong pos2kmer (tref_cat_binary_t *btref, int nllong, int strnd, int pos, int K) {
	llong temp = 0LLU;
	int ind = pos/32 ; //which element in cat_tref 	
	int rmd = (pos % 32)*2; //# reminder in bits
	int lshft = (32-K)*2; //#left empty bits

	if(strnd == 1){ //use positive strnd
		temp = ( (btref->cat_tref)[ind] << rmd ) >> lshft ;
		if( rmd > lshft ) 			
			temp += (btref->cat_tref)[ind+1] >> (64 - rmd + lshft); 
	}
	else{
		temp = ( (btref->rc_cat_tref)[nllong - 1 - ind] >>rmd) & (_64MASK >> lshft)  ;

		if(rmd > lshft )
			temp +=	 ((btref->rc_cat_tref)[nllong - 2 - ind] << (64-(rmd-lshft)) ) >> lshft;
	}
	
	return (temp);

}
*/
ext_btref_kmer_index_t*  build_ext_btref_kmer_index (btref_kmer_index_t* btref_kmer_index, tref_cat_binary_t *btref, int K){
	assert(K<=32);
	ext_btref_kmer_index_t* ext_btref_kmer_index = malloc(sizeof(ext_btref_kmer_index_t));
	int htsize = nextPrime( (int)((double)btref->n / LD_FCTR ) ); //size of hashtable
	printf("htsize=%d\n",htsize);
	int *kmer_pos_ht = malloc(htsize * sizeof(int)); //hashtable for kmer 1st present position
	for(int i=0;i<htsize;i++) kmer_pos_ht[i] = HASH_DFLT;	

	int kmer_ht_rank = 0 ;
	int **kmer_pos_mtx = NULL;
	llong kmer, crkmer,minkmer;	
	int nllong = btref->n % 32 == 0 ?  btref->n/32 :  (btref->n/32 + 1); //nelement

	int uniq_kmer_n = 0;
	int npos = 0;
	for(int p = 0; p < btref->n; p++){
		if(btref_kmer_index[p].kmer_hash_rank == POS_INVLD ) continue;
		npos++;
		kmer = pos2kmer(btref,nllong,1,p,K) ;
		crkmer =  pos2kmer(btref,nllong,-1,p,K) ;	
		int strnd = 1;// default postive strand

		if(crkmer < kmer) {
			minkmer = crkmer;
			strnd = -1 ;
		}
		else {
			minkmer = kmer;
			strnd = 1;
		}
		//hash
		for(int i = 0; i < htsize; i++){
			int hsv = HASH(minkmer,i,htsize);
			if(kmer_pos_ht[hsv] == HASH_DFLT){
				kmer_pos_ht[hsv] = strnd * (p+1);
				uniq_kmer_n++;
				break;
			}
			else{
				int cld_hs	= kmer_pos_ht[hsv] ;
				int cld_hs_strnd = cld_hs > 0? 1 : -1;
				int abs_fst_pos = abs(cld_hs) - 1;				
				llong cld_kmer = pos2kmer(btref,nllong,cld_hs_strnd,abs_fst_pos, K) ;
				if(cld_kmer  == minkmer){
					if( btref_kmer_index[abs_fst_pos].kmer_hash_rank == KMER_HSH_RNK_DFLT){
						btref_kmer_index[abs_fst_pos].kmer_hash_rank = kmer_ht_rank;//namely kmer_pos_mtx row index
		
#define REALC_BLK 1024 
						if(kmer_ht_rank % REALC_BLK == 0)	
							kmer_pos_mtx = realloc(kmer_pos_mtx,( kmer_ht_rank / REALC_BLK + 1) * REALC_BLK * sizeof(int *)) ;					
						
						kmer_pos_mtx[kmer_ht_rank] = malloc(2*sizeof(int ));
						kmer_pos_mtx[kmer_ht_rank][0] = 2; //array size
						kmer_pos_mtx[kmer_ht_rank][1] = strnd * (p+1);;
						kmer_ht_rank++;//namely kmer_pos_mtx row index
					}
					else{

						int *row = kmer_pos_mtx[btref_kmer_index[abs_fst_pos].kmer_hash_rank];
						row = realloc(row,(row[0]+1)*sizeof(int ));
						kmer_pos_mtx[btref_kmer_index[abs_fst_pos].kmer_hash_rank] = row;
						row[row[0]] = strnd * (p+1);
						row[0]++;	 	 
					}            										
					break;
				}
				else continue;
			}
		}
	};

	/* convert kmer_pos_ht to tref_uniq_kmer_pos array (reference kmer(pos) array),
	 which defined the rows of A and b (2022-5-27) */
	int *tref_uniq_kmer_pos = malloc(uniq_kmer_n * sizeof(int));
    int j = 0;
    for (int i = 0 ; i < htsize; i++){
        if( kmer_pos_ht[i] != HASH_DFLT){
            tref_uniq_kmer_pos[j] = kmer_pos_ht[i];
            j++;
        }
    }
  free(kmer_pos_ht);

	ext_btref_kmer_index->K = K;
	ext_btref_kmer_index->htsize = htsize;
	ext_btref_kmer_index->uniq_kmer_n = uniq_kmer_n;
	ext_btref_kmer_index->npos = npos; //# of valid pos
	ext_btref_kmer_index->btref_kmer_index = btref_kmer_index;
	ext_btref_kmer_index->ref = tref_uniq_kmer_pos; //kmer_pos_ht;(2022-5-27)
	ext_btref_kmer_index->kmer_pos_mtx =  kmer_pos_mtx;
	return(ext_btref_kmer_index);
}

static inline int t_pos2wght (int t_pos,int tlen, int fragl, int readl, int K) {
	//a = 0; b = tlen - fragl;
	int b = tlen - fragl;
	int c = t_pos + K - readl;
    //int d = t_pos;
    int cp = t_pos + K - fragl;
    int dp = t_pos + readl - fragl;
	
    int A = 0 < c ? c : 0 ;
    int B = b< t_pos ? b: t_pos ;
    int Ap = 0 < cp? cp:0;
    int Bp = b < dp ? b : dp;

	int w1 = B - A + 1;
	if(w1<0) w1 = 0 ;

    int w2 = Bp - Ap + 1  ;
    if(w2<0) w2 = 0 ;
	
	return (w1 + w2);
}

/*build_Kref_mtx and free  ext_btref_kmer_index internal*/
Kref_mtx_t* build_Kref_mtx (ext_btref_kmer_index_t*  ext_btref_kmer_index, tref_cat_binary_t *btref, indexProperty_t * indexProperty){

    int K = indexProperty->K;
//    int htsize = ext_btref_kmer_index->htsize;
    int uniq_kmer_n = ext_btref_kmer_index->uniq_kmer_n;
    int npos = ext_btref_kmer_index->npos;
    int noT = btref->noT;
    btref_kmer_index_t * btref_kmer_index = ext_btref_kmer_index->btref_kmer_index;
    // tref_uniq_kmer_pos (ref) defines the rows of A and b*/
		int *tref_uniq_kmer_pos = ext_btref_kmer_index->ref;
    int **kmer_pos_mtx = ext_btref_kmer_index->kmer_pos_mtx;

    int *Tsum_wghts = calloc(noT, sizeof(int));
    int *kmer_tids = malloc(npos*sizeof(int));
    double *kmer_wgts = malloc(npos*sizeof(double));
    int *kmer_rowp = malloc( (uniq_kmer_n + 1) * sizeof(int)) ;
    kmer_rowp[0] = 0; // int kmer_rowp_ind = 0;

    int n_wgts = 0; //current row, current column (tid)'s weight, check if equal to # of non-zero
    int fragl = indexProperty->fragl;
    int readl = indexProperty->readl;
    int abs_pos,tid,t_start,tlen,t_pos,wght;

    for (int i = 0; i < uniq_kmer_n; i++ )  {
        abs_pos = abs(tref_uniq_kmer_pos[i]) - 1 ;
        tid = btref_kmer_index[abs_pos].tid;
        t_start = tid == 0 ? 0 : (btref->tref_nxt_start)[tid - 1];
        tlen = (btref->tref_nxt_start)[tid] - t_start;
        t_pos = abs_pos - t_start; //position to this transcript start
        wght = t_pos2wght(t_pos,tlen,fragl,readl,K);
/*initilize for each uniq_kmer, either rank == KMER_HSH_RNK_DFLT or not*/
        Tsum_wghts[tid] += wght;
        kmer_tids[n_wgts] = tid;
        kmer_wgts[n_wgts] = wght;
        kmer_rowp[i+1] = kmer_rowp[i] + 1; //i+1th row start pointer

        int rank = btref_kmer_index[abs_pos].kmer_hash_rank ;
        if(rank != KMER_HSH_RNK_DFLT) { //duplicated kmer
            int *row = kmer_pos_mtx[rank];
            // get pos in this row, row[0] is arrsize, shoule not be counted
            for(int p = 1; p < row[0]; p++){
                int p_abs_pos = abs(row[p]) - 1 ;
                int p_tid = btref_kmer_index[p_abs_pos].tid;
                int p_t_start = p_tid == 0 ? 0 : (btref->tref_nxt_start)[p_tid - 1];
                int p_tlen = (btref->tref_nxt_start)[p_tid] - p_t_start;
                int p_t_pos = p_abs_pos - p_t_start;
                int p_wght = t_pos2wght(p_t_pos,p_tlen,fragl,readl,K);

                if(p_tid != tid){ // current tid ended
                    kmer_rowp[i+1]++;// ith kmer row found one more tid, so i+1 th row start increase 1;
                    n_wgts++;//begin next column(tid), so increase n_wgts
                    kmer_tids[n_wgts] = p_tid;
                    tid = p_tid;
                }

                kmer_wgts[n_wgts] += p_wght;
                Tsum_wghts[p_tid] += p_wght;
            }
			free(row); //free kmer_pos_mtx:, make sure will not access this row again
        };
        n_wgts++; // begin next row, so increase n_wgts
    }
	free(btref_kmer_index);
	free(kmer_pos_mtx); // all ext_btref_kmer_index freed
	printf("\t< kmer_Twght_mtx builded, free btref_kmer_index and kmer_pos_mtx >\t%dM released\n",btref->n * 8/1024/1024);
    /*normalize kmer_wgts[i] by ith trnscript total wght*/

    for(int i = 0; i < n_wgts ; i++){
        int temp_tid = kmer_tids[i];
		int temp_t_stat = temp_tid == 0? 0:(btref->tref_nxt_start)[temp_tid - 1];
        int temp_tlen = (btref->tref_nxt_start)[temp_tid] - temp_t_stat;
		//caution: non ACGT also included here, may need modify
        //version 0: kmer_wgts[i] = (double)kmer_wgts[i] / Tsum_wghts[temp_tid]* (temp_tlen - K + 1) ; 
		//version 1: kmer_wgts[i] = (double)kmer_wgts[i] / Tsum_wghts[temp_tid] * (temp_tlen - K + 1) * (temp_tlen - fragl + 1) / temp_tlen ; 
		//version 2: 
		int mod = temp_tlen % fragl;  
		kmer_wgts[i] = (double)kmer_wgts[i] / (double)Tsum_wghts[temp_tid] * (double) (temp_tlen - mod - (1 +  temp_tlen / fragl)*K) ;
    };

	free(Tsum_wghts);

	kmer_Twght_mtx_t* kmer_Twght_mtx = malloc(sizeof(kmer_Twght_mtx_t));
	kmer_Twght_mtx->nzmax = n_wgts;
	kmer_Twght_mtx->nrow = uniq_kmer_n;
	kmer_Twght_mtx->ncol = noT;
	kmer_Twght_mtx->p = kmer_rowp;
	kmer_Twght_mtx->i = kmer_tids;
	kmer_Twght_mtx->x = kmer_wgts;
	//definition of the row of kmer_Twght_mtx 
	Kref_t* Kref = malloc(sizeof(Kref_t)); 	
	Kref->K = K;	
	Kref->size = uniq_kmer_n;
	Kref->ref = tref_uniq_kmer_pos;

	Kref_mtx_t *Kref_mtx = malloc(sizeof(Kref_mtx_t));//kmer_Twght_mtx +  Kref->ref
	Kref_mtx->kmer_Twght_mtx = kmer_Twght_mtx;
	Kref_mtx->Kref = Kref;	
			
	printf("\t< ReferenceKmer-TransctriptID sparse Matrix builded >\tsize: %d x %d, nnz=%d\n",uniq_kmer_n,noT,n_wgts );
	return(Kref_mtx);
}
/*convert kmer_Twght_mtx to cholmod_sparse compatible matrix, row to column; */
cholmod_sparse* kmer_Twght_mtx2cholmod_sparse(kmer_Twght_mtx_t *kmer_Twght_mtx){

	cholmod_sparse* M = malloc(sizeof(cholmod_sparse));
	M->nzmax = kmer_Twght_mtx->nzmax ;
	M->nrow = kmer_Twght_mtx->ncol;
	M->ncol = kmer_Twght_mtx->nrow;
	M->p = kmer_Twght_mtx->p;
	M->i = kmer_Twght_mtx->i;
	M->x = kmer_Twght_mtx->x;
	M->stype = 0; // matrix is "unsymmetric"
	M->itype = CHOLMOD_INT;
	M->xtype = CHOLMOD_REAL; // pattern, real, complex, or zomplex
	M->dtype = 1;// x and z are double or float
	M->sorted = 0;// TRUE if columns are sorted, FALSE otherwise
	M->packed = 1;// TRUE if packed (nz ignored), FALSE if unpacked (nz is required)
	return(M);
}

pos2bidx_t *build_pos2bidx(ext_btref_kmer_index_t*  ext_btref_kmer_index, tref_cat_binary_t *btref){
	
	pos2bidx_t *Pos2bidx = malloc(sizeof(pos2bidx_t));	
	// set array pos2bidx default value to 0, may need use other value 
	int* pos2bidx = calloc( btref->n, sizeof(int) );
	
	int uniq_kmer_n = ext_btref_kmer_index->uniq_kmer_n;

	btref_kmer_index_t * btref_kmer_index = ext_btref_kmer_index->btref_kmer_index;
	int *tref_uniq_kmer_pos = ext_btref_kmer_index->ref;
	int **kmer_pos_mtx = ext_btref_kmer_index->kmer_pos_mtx;
	
	for (int i = 0; i < uniq_kmer_n; i++ )  {
		int abs_pos = abs(tref_uniq_kmer_pos[i]) - 1 ;
		pos2bidx[abs_pos]	= i;
		
		int rank = btref_kmer_index[abs_pos].kmer_hash_rank ;	
		if(rank != KMER_HSH_RNK_DFLT) {
			int *row = kmer_pos_mtx[rank];
			// get pos in this row, row[0] is arrsize, shoule not be counted						
      for(int p = 1; p < row[0]; p++)
				pos2bidx[abs(row[p]) - 1] = i ;						
		}	
	}
	Pos2bidx->len = btref->n;
	Pos2bidx->pos2bidx = pos2bidx;
	return Pos2bidx;
}



