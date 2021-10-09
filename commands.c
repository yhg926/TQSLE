#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "basic.h"
#include "index_core.h"
#include "quant_core.h"
#include "IO.h"
#include "cholmod.h"
#include "cholmod_core.h"


char path_holder [256];

int index_cmd (indexProperty_t* indexProperty ){

    struct timeval t1, t2;
    double elapsedTime;
	
	assert ( indexProperty->fragl > indexProperty->readl ); assert ( indexProperty->readl > indexProperty->K );
	assert(indexProperty->tref_fa != NULL);
	// references shorter than indexProperty->fragl will be flitered 
	trscrpt_ref_arr_t *tref = Read_fas_trscrpt(indexProperty->tref_fa, indexProperty->fragl);
		
	sprintf(path_holder,"%s/%s",indexProperty->outpath, tnamelen_f);
	fprint_tref_namelen (path_holder, tref);
	
	tref_cat_binary_t * btref = cat_tref_arr(tref) ;
	trscrpt_ref_arr_free(tref);
	
	btref_kmer_index_t* btref_kmer_index = build_btref_kmer_index ( btref,  indexProperty->K);
	printf("btref_kmer_index builded\n");

	gettimeofday(&t1, NULL);
	ext_btref_kmer_index_t * ext_btref_kmer_index = build_ext_btref_kmer_index (btref_kmer_index, btref,  indexProperty->K);
	gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("ext_btref_kmer_index builded. Time: %.4f ms\n",elapsedTime);
	
	Kref_mtx_t* Kref_mtx = build_Kref_mtx(ext_btref_kmer_index, btref, indexProperty);
	gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("Kref_mtx builded and ext_btref_kmer_index freed. Time: %.4f ms\n",elapsedTime);

	sprintf(path_holder,"%s/%s",indexProperty->outpath, btref_f); write_btref(path_holder,btref); free_btref(btref);
	sprintf(path_holder,"%s/%s",indexProperty->outpath, Kref_f); write_Kref(path_holder, Kref_mtx->Kref);free(Kref_mtx->Kref);
	sprintf(path_holder,"%s/%s",indexProperty->outpath, kmer_Twght_mtx_f); write_kmer_Twght_mtx(path_holder,Kref_mtx->kmer_Twght_mtx);
	gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("btref, Kref and kmer_Twght_mtx written, converting to At. Time: %.4f ms\n",elapsedTime);
	
	cholmod_sparse* At = kmer_Twght_mtx2cholmod_sparse(Kref_mtx->kmer_Twght_mtx);
	
	
    cholmod_common* c = (cholmod_common*)malloc(sizeof(cholmod_common));
    cholmod_start (c) ;

	cholmod_sparse* AtA = cholmod_aat(At,NULL,0,1,c);	
	
	gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("AtA builed. Time: %.4f ms\n",elapsedTime);
	free(At->x);free(At->i);free(At->p);free(At);	
	cholmod_sparse* I = cholmod_speye(AtA->nrow,AtA->ncol,AtA->xtype,c);
	double alpha[2] = {1,0};
	double beta[2] = {10,0};
	cholmod_sparse* AtAI = cholmod_add(AtA,I,alpha,beta,1,0,c);
	AtAI->stype = 1;
	cholmod_free_sparse(&AtA,c);
	gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("AtAI created.Time: %.4f ms\n",elapsedTime);	

	cholmod_factor *L = cholmod_analyze (AtAI, c) ;
    gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;	
	printf("cholmod_analyzed,Time: %.4f ms\n",elapsedTime);
	cholmod_factorize (AtAI, L, c) ;
	gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("cholmod_factorized,Time: %.4f ms\n",elapsedTime);
	
	sprintf(path_holder,"%s/%s",indexProperty->outpath, factorL_f); cholmod_factor_write(path_holder,L);
    gettimeofday(&t2, NULL); elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
	printf("cholmod_factor writed,Time: %.4f ms\n",elapsedTime);
	
	return(0);
}

int quant_cmd (quantProperty_t * quantProperty) {


	sprintf(path_holder,"%s/%s",quantProperty->indexpath, btref_f);
    tref_cat_binary_t * btref = read_btref(path_holder);
    printf("noT=%d\n",btref->noT);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, Kref_f);
    Kref_t *Kref = read_Kref(path_holder);

    basis_Kref_t* basis_Kref = build_basisKref(btref,Kref);
    printf("basis_Kref loaded\n");

    double *db = create_b (btref, basis_Kref, quantProperty->num_remaining_args, quantProperty->remaining_args);

    cholmod_common *c = (cholmod_common*)malloc(sizeof(cholmod_common));
    cholmod_start(c);
    int size = basis_Kref->Kref_size ;
    cholmod_dense *b = make_cholmod_dense(db, size , c);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, kmer_Twght_mtx_f);
    kmer_Twght_mtx_t* kmer_Twght_mtx = read_kmer_Twght_mtx(path_holder);
    cholmod_sparse* At = kmer_Twght_mtx2cholmod_sparse(kmer_Twght_mtx);

    cholmod_dense *Atb = cholmod_zeros(At->nrow,1,CHOLMOD_REAL,c);

    double alpha[2] = {1,0};
    double beta[2] = {0,0};
    cholmod_sdmult(At, 0, alpha ,  beta , b , Atb, c);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, factorL_f);
    cholmod_factor *L = cholmod_factor_read(path_holder,c);
    printf("L loaded \n");
    cholmod_dense *rslt = cholmod_solve (CHOLMOD_A, L, Atb, c) ;
    printf("x get sloved \n");
	cholmod_finish(c);	

	//printf results
	fprint_abundance (quantProperty, rslt->nrow, rslt->x);

    return(1);
}



















