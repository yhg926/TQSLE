#ifndef IO_H
#define IO_H

#include "index_core.h"

extern const char tnamelen_f[];
extern const char btref_f[] ;
extern const char Kref_f[] ;
extern const char kmer_Twght_mtx_f[] ;
extern const char factorL_f[] ;
extern const char abundance_f[] ;
extern const char pos2bidx_f[] ; //btref sequence position to b index 

int write_pos2bidx ( char* , pos2bidx_t *);
pos2bidx_t *read_pos2bidx (char*);
int fprint_tref_namelen (char *, trscrpt_ref_arr_t *);
int fprint_abundance (quantProperty_t *, int, double *, double);
int write_btref( char *, tref_cat_binary_t * );
tref_cat_binary_t * read_btref(char *) ;

int write_Kref(char*, Kref_t* );
Kref_t * read_Kref(char*);

int write_kmer_Twght_mtx(char*, kmer_Twght_mtx_t*);
kmer_Twght_mtx_t* read_kmer_Twght_mtx(char *);

int cholmod_factor_write(char*, cholmod_factor *);
cholmod_factor *cholmod_factor_read(char*, cholmod_common *);

#endif





