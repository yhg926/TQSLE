#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "basic.h"
#include "cholmod.h"
#include "cholmod_internal.h"
#include "cholmod_core.h"
#include "index_core.h"
#include "IO.h"

const char tnamelen_f[]= "tnamelen.out";
const char btref_f[] = "btref.out";
const char Kref_f[] = "Kref.out";
const char kmer_Twght_mtx_f[] = "kmer_Twght_mtx.out";
const char factorL_f[] = "factorL.out";
const char abundance_f[] = "abundance.out";

int fprint_tref_namelen (char *tnamelen_fout, trscrpt_ref_arr_t *tref_arr){
	FILE *fh = fopen(tnamelen_fout, "w");
	assert( (fh != NULL) && "File open failed in write_tref_namelen ()" );

	for(int i = 0; i < tref_arr->n ; i++)
		fprintf(fh,"%s\t%d\n",(tref_arr->trscrpt_ref_idx)[i].name,(tref_arr->trscrpt_ref_idx)[i].len);

	return(1);
}

int fprint_abundance (quantProperty_t *quantProperty, int nrow, double *x, double neg_fix){

    char line [1024];
	int num_neg_fix = 0 ;
	int x_sum  = 0;
	for(int i = 0; i< nrow ; i++ ) x_sum += x[i] ;

    sprintf(line,"%s/%s",quantProperty->outpath, abundance_f);
    FILE* fout = fopen(line, "w"); assert(fout!=NULL);

    sprintf(line,"%s/%s",quantProperty->indexpath, tnamelen_f);
    FILE* tname_fh = fopen( line  , "r"); assert(tname_fh!=NULL);

    fprintf(fout,"ID\tLength\tAbundance\tTPM\n");

    for(int i = 0; i< nrow ; i++ ){
        fgets(line, 1024, tname_fh);
		*(strchr(line, '\n')) = '\0';
		
		if( neg_fix != 0 ) {

			if (x[i] < 0){
				x[i] = neg_fix == 1 ? 0 : neg_fix ; 
				num_neg_fix++;
			}
		}

        fprintf(fout,"%s\t%lf\t%lf\n",line, x[i], (double)x[i] / x_sum * 1000000 );
    }

    fclose(fout);
    fclose(tname_fh);
	return(num_neg_fix);
}


//---------binary tref seq----------
int write_btref( char *btref_fout, tref_cat_binary_t * btref)
{
    FILE *fh = fopen(btref_fout, "wb");
	assert( (fh != NULL) && "File open failed in write_btref ()" ) ; 

    fwrite(&(btref->n),sizeof(int),1,fh);
    fwrite(btref->cat_tref,8, (int)(btref->n /32) + 1 , fh);

    fwrite(&(btref->noN), sizeof(int), 1, fh);
    fwrite(btref->Npos, sizeof(int), btref->noN, fh);

    fwrite(&(btref->noT), sizeof(int), 1, fh);
    fwrite(btref->tref_nxt_start, sizeof(int), btref->noT, fh);

    fclose(fh);
    return(1);
}

tref_cat_binary_t * read_btref( char *btref_fin){
    FILE *fh = fopen(btref_fin, "rb");
    assert( (fh != NULL) && "File open failed in read_btref ()" ) ;

    tref_cat_binary_t * btref = malloc(sizeof(tref_cat_binary_t));

    fread( &(btref->n) ,sizeof(int),1,fh);
    btref->cat_tref =  (llong *)malloc( ((int)(btref->n /32) + 1)*sizeof(llong) ) ;
    fread(btref->cat_tref,8, (int)(btref->n /32) + 1 , fh);

    fread(&(btref->noN), sizeof(int), 1, fh);
    btref->Npos = (int *)malloc(btref->noN * sizeof(int) );
    fread(btref->Npos, sizeof(int), btref->noN, fh);


    fread(&(btref->noT), sizeof(int), 1, fh);
    btref->tref_nxt_start = (int *)malloc(btref->noT * sizeof(int) );
    fread(btref->tref_nxt_start, sizeof(int), btref->noT, fh);

    btref->rc_cat_tref = NULL;
    crvs_btref(btref) ;

    fclose(fh);
    return(btref);
}
//-------------ref Kmer: K + ref pos(definition of kmer_Twght_mtx rows)-----------------
int write_Kref( char *Kref_f, Kref_t* Kref)
{
    FILE *fh = fopen(Kref_f, "wb");
    assert( (fh != NULL) && "File open failed in write_Kref()") ;

    fwrite(&(Kref->K),sizeof(int),1,fh);
    fwrite(&(Kref->size),sizeof(int),1,fh);
    fwrite(Kref->ref,sizeof(int),Kref->size,fh);
    fclose(fh);
    return(1);
}

Kref_t * read_Kref(char *Kref_f){
    FILE *fh = fopen(Kref_f, "rb");
    assert( (fh != NULL) && "File open failed in read_Kref()") ;
    Kref_t *Kref = malloc(sizeof(Kref_t));
    fread(&(Kref->K), sizeof(int) , 1, fh);
    fread(&(Kref->size),sizeof(int),1,fh);
    Kref->ref = malloc( Kref->size * sizeof(int));
    fread(Kref->ref,sizeof(int),Kref->size,fh);
    fclose(fh);
    return(Kref);
}
//-----------------Matrix A: kmer_Twght_mtx-------
int write_kmer_Twght_mtx( char *kmer_Twght_mtx_f, kmer_Twght_mtx_t *kmer_Twght_mtx){
    FILE *fh = fopen(kmer_Twght_mtx_f, "wb");
    assert((fh != NULL) && "File open failed in write_kmer_Twght_mtx()") ;
    fwrite(&(kmer_Twght_mtx->nzmax), sizeof(int) ,1, fh);
    fwrite(&(kmer_Twght_mtx->nrow), sizeof(int) ,1, fh);
    fwrite(&(kmer_Twght_mtx->ncol), sizeof(int) ,1, fh);
    fwrite(kmer_Twght_mtx->p, sizeof(int), kmer_Twght_mtx->nrow + 1, fh);
    fwrite(kmer_Twght_mtx->i, sizeof(int), kmer_Twght_mtx->nzmax, fh);
    fwrite(kmer_Twght_mtx->x, sizeof(double), kmer_Twght_mtx->nzmax, fh);
    fclose(fh);
    return(1);
}

kmer_Twght_mtx_t* read_kmer_Twght_mtx(char *kmer_Twght_mtx_f){

    FILE *fh = fopen(kmer_Twght_mtx_f, "rb");
    assert((fh != NULL)&& "File open failed in read_kmer_Twght_mtx()") ;
    kmer_Twght_mtx_t* kmer_Twght_mtx = malloc(sizeof(kmer_Twght_mtx_t));
    fread(&(kmer_Twght_mtx->nzmax), sizeof(int) ,1, fh);
    fread(&(kmer_Twght_mtx->nrow), sizeof(int) ,1, fh);
    fread(&(kmer_Twght_mtx->ncol), sizeof(int) ,1, fh);

    kmer_Twght_mtx->p = malloc(sizeof(int)*(kmer_Twght_mtx->nrow + 1));
    fread(kmer_Twght_mtx->p, sizeof(int), kmer_Twght_mtx->nrow + 1, fh);

    kmer_Twght_mtx->i = malloc( sizeof(int) * kmer_Twght_mtx->nzmax );
    fread(kmer_Twght_mtx->i, sizeof(int), kmer_Twght_mtx->nzmax, fh);

    kmer_Twght_mtx->x = malloc(sizeof(double) * kmer_Twght_mtx->nzmax);
    fread(kmer_Twght_mtx->x, sizeof(double), kmer_Twght_mtx->nzmax, fh);

    fclose(fh);
    return(kmer_Twght_mtx);
}

/*cholmod factor w/r*/
int CHOLMOD(factor_write)
(
     char *factorL_f, /* file to write to, must already be open */
    cholmod_factor *L
)
{
	FILE *f = fopen(factorL_f, "wb");
    assert((f != NULL) && "File  open failed in CHOLMOD(factor_write)()") ;	

    /*Write all non-pointer elements */
    fwrite( &(L->n), sizeof(L->n), 1, f); //1

    fwrite( &(L->minor), sizeof(L->minor), 1, f); //2
    fwrite( &(L->nzmax), sizeof(L->nzmax), 1, f); //3
    fwrite( &(L->nsuper), sizeof(L->nsuper), 1, f); //4
    fwrite( &(L->ssize), sizeof(L->ssize), 1, f); //5
    fwrite( &(L->xsize), sizeof(L->xsize), 1, f); //6
    fwrite( &(L->maxcsize), sizeof(L->maxcsize), 1, f); //7
    fwrite( &(L->maxesize), sizeof(L->maxesize), 1, f); //8
    fwrite( &(L->ordering), sizeof(L->ordering), 1, f);//9
    fwrite( &(L->is_ll), sizeof(L->is_ll), 1, f); //10
    fwrite( &(L->is_super), sizeof(L->is_super), 1, f); //11
    fwrite( &(L->is_monotonic), sizeof(L->is_monotonic), 1, f); //12
    fwrite( &(L->itype), sizeof(L->itype), 1, f); //13
    fwrite( &(L->xtype), sizeof(L->xtype), 1, f); //14
    fwrite( &(L->dtype), sizeof(L->dtype), 1, f); //15
    fwrite( &(L->useGPU), sizeof(L->useGPU), 1, f); //16

    size_t n =  L->n;
    size_t nzmax = L->nzmax;
    size_t nsuper = L->nsuper ;

    /*Write all pointer elements*/
    fwrite( L->Perm, sizeof(Int), n, f);//17
    fwrite( L->ColCount, sizeof(Int), n, f);//18


    if (!(L->is_super)) //origin: if (L->xtype != CHOLMOD_PATTERN && !(L->super))
    {
    fwrite( L->p, sizeof(Int), n+1, f);//19
    fwrite( L->prev, sizeof(Int), n+2, f);//20
    fwrite( L->next, sizeof(Int), n+2, f);//21
    fwrite( L->nz, sizeof(Int), n, f);//22
    fwrite( L->i, sizeof(Int), nzmax, f);//23
    fwrite( L->x, sizeof(double), nzmax, f);//24
    }
    else if(L->is_super)
    {
    fwrite( L->super, sizeof(Int), nsuper+1, f);//25
    fwrite( L->pi, sizeof(Int), nsuper+1, f);//26
    fwrite( L->px, sizeof(Int), nsuper+1, f);//27
    fwrite( L->s, sizeof(Int), L->ssize, f);//28
    fwrite( L->x, sizeof(double), L->xsize, f);//29
    }

	fclose(f);

    return 1;
};

cholmod_factor *CHOLMOD(factor_read)
(
    char *factorL_f, /* file to read, must already be open */
    cholmod_common *Common
)
{
	FILE *f = fopen(factorL_f, "rb");
	assert( (f != NULL) && "File  open failed in CHOLMOD(factor_read)()")  ;
    size_t n = 0 ;
    assert(fread( &n, sizeof(n), 1, f) == 1) ;//1
    cholmod_factor *L = CHOLMOD(allocate_factor) (n, Common) ;
    assert(sizeof(n) == sizeof(L->n));//make sure L->n type is same size to n
    L->n = n;
    /*Read all non-pointer elements */
    assert(fread( &(L->minor), sizeof(L->minor), 1, f) == 1); //2
    assert(fread( &(L->nzmax), sizeof(L->nzmax), 1, f) == 1); //3
    assert(fread( &(L->nsuper), sizeof(L->nsuper), 1, f) == 1); //4
    assert(fread( &(L->ssize), sizeof(L->ssize), 1, f) == 1); //5
    assert(fread( &(L->xsize), sizeof(L->xsize), 1, f) == 1); //6
    assert(fread( &(L->maxcsize), sizeof(L->maxcsize), 1, f) == 1); //7
    assert(fread( &(L->maxesize), sizeof(L->maxesize), 1, f) == 1); //8
    assert(fread( &(L->ordering), sizeof(L->ordering), 1, f) == 1);//9
    assert(fread( &(L->is_ll), sizeof(L->is_ll), 1, f) == 1); //10
    assert(fread( &(L->is_super), sizeof(L->is_super), 1, f) == 1); //11
    assert(fread( &(L->is_monotonic), sizeof(L->is_monotonic), 1, f) == 1); //12
    assert(fread( &(L->itype), sizeof(L->itype), 1, f) == 1); //13
    assert(fread( &(L->xtype), sizeof(L->xtype), 1, f) == 1); //14
    assert(fread( &(L->dtype), sizeof(L->dtype), 1, f) == 1); //15
    assert(fread( &(L->useGPU), sizeof(L->useGPU), 1, f) == 1); //16

    size_t nzmax = L->nzmax;
    size_t nsuper = L->nsuper ;

    assert(fread( L->Perm, sizeof(Int), n, f) == n);//17
    assert(fread( L->ColCount, sizeof(Int), n, f) == n);//18


    if (!(L->is_super)) //origin: if (L->xtype != CHOLMOD_PATTERN && !(L->super))
    {
    L->p = malloc(sizeof(Int)*(n+1));
    assert(fread( L->p, sizeof(Int), n+1, f) == (n+1) ) ;//19

    L->prev = malloc(sizeof(Int)*(n+2));
    assert(fread( L->prev, sizeof(Int), n+2, f) == (n+2) );//20

    L->next = malloc(sizeof(Int)*(n+2));
    assert(fread( L->next, sizeof(Int), n+2, f) == (n+2) );//21

    L->nz = malloc(sizeof(Int)*n);
    assert(fread( L->nz, sizeof(Int), n, f) == n);//22

    L->i = malloc(sizeof(Int)*nzmax);
    assert(fread( L->i, sizeof(Int), nzmax, f) == nzmax);//23

    L->x = malloc(sizeof(double)*nzmax);
    assert(fread( L->x, sizeof(double), nzmax, f) == nzmax);//24
    }

    else if(L->is_super)
    {

    L->super = malloc(sizeof(Int)*(nsuper+1));
    assert(fread( L->super, sizeof(Int), nsuper+1, f) == (nsuper+1) );//25

    L->pi = malloc(sizeof(Int)*(nsuper+1));
    assert(fread( L->pi, sizeof(Int), nsuper+1, f) == (nsuper+1) );//26

    L->px = malloc(sizeof(Int)*(nsuper+1));
    assert(fread( L->px, sizeof(Int), nsuper+1, f) == (nsuper+1) );//27

    L->s = malloc(sizeof(Int) * (L->ssize+1));
    assert(fread( L->s, sizeof(Int), L->ssize, f) == L->ssize);//28

    L->x = malloc(sizeof(double) * (L->xsize+1));
    assert(fread( L->x, sizeof(double), L->xsize, f) == L->xsize);//29

    }
	fclose(f);
    return (L) ;
};











