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
    struct timeval t1, t2; gettimeofday(&t1, NULL);
	assert ( indexProperty->fragl >= indexProperty->readl ); assert ( indexProperty->readl >= indexProperty->K );
	assert(indexProperty->tref_fa != NULL);

	// references shorter than indexProperty->fragl will be flitered 
	trscrpt_ref_arr_t *tref = Read_fas_trscrpt(indexProperty->tref_fa, indexProperty->fragl);
		
	sprintf(path_holder,"%s/%s",indexProperty->outpath, tnamelen_f);
	fprint_tref_namelen (path_holder, tref);
	
	tref_cat_binary_t * btref = cat_tref_arr(tref) ;
	trscrpt_ref_arr_free(tref);
	
	btref_kmer_index_t* btref_kmer_index = build_btref_kmer_index ( btref,  indexProperty->K);
	gettimeofday(&t2, NULL); printf("\t< btref_kmer_index builded >\t%ld s\n", t2.tv_sec - t1.tv_sec);
	
	ext_btref_kmer_index_t * ext_btref_kmer_index = build_ext_btref_kmer_index (btref_kmer_index, btref,  indexProperty->K);
	gettimeofday(&t2, NULL); printf("\t< ext_btref_kmer_index builded >\t%ld s\n", t2.tv_sec - t1.tv_sec);
	
	pos2bidx_t *Pos2bidx =	build_pos2bidx(ext_btref_kmer_index,btref);
	sprintf(path_holder,"%s/%s",indexProperty->outpath, pos2bidx_f);
	write_pos2bidx(path_holder,Pos2bidx);
	free(Pos2bidx->pos2bidx); free(Pos2bidx) ;
	gettimeofday(&t2, NULL); printf("\t< Pos2bidx written >\t%ld s\n", t2.tv_sec - t1.tv_sec);
	

	Kref_mtx_t* Kref_mtx = build_Kref_mtx(ext_btref_kmer_index, btref, indexProperty);
	gettimeofday(&t2, NULL); printf("\t< Kref_mtx builded and ext_btref_kmer_index freed >\t%ld s\n", t2.tv_sec - t1.tv_sec);

	sprintf(path_holder,"%s/%s",indexProperty->outpath, btref_f); write_btref(path_holder,btref); free_btref(btref);
	sprintf(path_holder,"%s/%s",indexProperty->outpath, Kref_f); write_Kref(path_holder, Kref_mtx->Kref);free(Kref_mtx->Kref);
	sprintf(path_holder,"%s/%s",indexProperty->outpath, kmer_Twght_mtx_f); write_kmer_Twght_mtx(path_holder,Kref_mtx->kmer_Twght_mtx);
	gettimeofday(&t2, NULL); printf("\t< btref, Kref and kmer_Twght_mtx written, converting to At >\t%ld s\n",t2.tv_sec - t1.tv_sec);
	
	cholmod_sparse* At = kmer_Twght_mtx2cholmod_sparse(Kref_mtx->kmer_Twght_mtx);
	
    cholmod_common* c = (cholmod_common*)malloc(sizeof(cholmod_common));
    cholmod_start (c) ;

	cholmod_sparse* AtA = cholmod_aat(At,NULL,0,1,c);	
	gettimeofday(&t2, NULL); printf("\t< AtA builed >\t%ld s\n", t2.tv_sec - t1.tv_sec);

	free(At->x);free(At->i);free(At->p);free(At);	
	cholmod_sparse* I = cholmod_speye(AtA->nrow,AtA->ncol,AtA->xtype,c);
	double alpha[2] = {1,0};
	double beta[2] = {indexProperty->lamda,0};
	cholmod_sparse* AtAI = cholmod_add(AtA,I,alpha,beta,1,0,c);
	AtAI->stype = 1;

	gettimeofday(&t2, NULL); printf("\t< AtAI created >\t%ld s\n", t2.tv_sec - t1.tv_sec);	
	cholmod_free_sparse(&AtA,c);
	
	cholmod_factor *L = cholmod_analyze (AtAI, c) ;
	gettimeofday(&t2, NULL); printf("\t< cholmod_analyzed >\t%ld s\n", t2.tv_sec - t1.tv_sec);
	cholmod_factorize (AtAI, L, c) ;
	gettimeofday(&t2, NULL); printf("\t< cholmod_factorized >\t%ld s\n", t2.tv_sec - t1.tv_sec);
	
	sprintf(path_holder,"%s/%s",indexProperty->outpath, factorL_f); cholmod_factor_write(path_holder,L);
	gettimeofday(&t2, NULL); printf("\t< cholmod_factor written >\t%ld s\n", t2.tv_sec - t1.tv_sec);
	
	return(0);
}

int quant_cmd (quantProperty_t * quantProperty) {
	struct timeval t1, t2; gettimeofday(&t1, NULL);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, btref_f);
    tref_cat_binary_t * btref = read_btref(path_holder);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, Kref_f);
    Kref_t *Kref = read_Kref(path_holder);

    basis_Kref_t* basis_Kref = build_basisKref(btref,Kref);
    gettimeofday(&t2, NULL); printf("\t< basis_Kref loaded >\t%ld s\n", t2.tv_sec - t1.tv_sec);

    //double *db = create_b (btref, basis_Kref, quantProperty->num_remaining_args, quantProperty->is_fa, quantProperty->remaining_args);

		sprintf(path_holder,"%s/%s",quantProperty->indexpath, pos2bidx_f);
		pos2bidx_t * Pos2bidx = read_pos2bidx(path_holder);
		gettimeofday(&t2, NULL); printf("\t< Pos2bidx readed >\t%ld s\n", t2.tv_sec - t1.tv_sec);
		double *db = create_b2 (Pos2bidx, btref, basis_Kref, quantProperty->num_remaining_args, quantProperty->is_fa, quantProperty->remaining_args);		
		gettimeofday(&t2, NULL); printf("\t< raw b created >\t%ld s\n", t2.tv_sec - t1.tv_sec);
		free(Pos2bidx->pos2bidx); free(Pos2bidx);

    cholmod_common *c = (cholmod_common*)malloc(sizeof(cholmod_common));
    cholmod_start(c);
    int size = basis_Kref->Kref_size ;
    cholmod_dense *b = make_cholmod_dense(db, size , c);
	gettimeofday(&t2, NULL); printf("\t< b created >\t%ld s\n", t2.tv_sec - t1.tv_sec);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, kmer_Twght_mtx_f);
    kmer_Twght_mtx_t* kmer_Twght_mtx = read_kmer_Twght_mtx(path_holder);
    cholmod_sparse* At = kmer_Twght_mtx2cholmod_sparse(kmer_Twght_mtx);

    cholmod_dense *Atb = cholmod_zeros(At->nrow,1,CHOLMOD_REAL,c);

    double alpha[2] = {1,0};
    double beta[2] = {0,0};
    cholmod_sdmult(At, 0, alpha ,  beta , b , Atb, c);
	gettimeofday(&t2, NULL); printf("\t< got Atb >\t%ld s\n", t2.tv_sec - t1.tv_sec);

	sprintf(path_holder,"%s/%s",quantProperty->indexpath, factorL_f);
    cholmod_factor *L = cholmod_factor_read(path_holder,c);
    printf("\t< L loaded >\n");
    cholmod_dense *rslt = cholmod_solve (CHOLMOD_A, L, Atb, c) ;
    gettimeofday(&t2, NULL); printf("\t< x get sloved >\t%ld s\n",  t2.tv_sec - t1.tv_sec);
	cholmod_finish(c);	

	//printf results (set neg_fix to 1 to correct negative value)
	int num_fix_neg = fprint_abundance (quantProperty, rslt->nrow, rslt->x, quantProperty->fix_neg);
		
	if (quantProperty->fix_neg != 0) {
		double xv = quantProperty->fix_neg == 1 ? 0 : quantProperty->fix_neg ;		
	printf("x: %d negative values were set to %lf\n", num_fix_neg, xv );
	}
    return(1);
}



















