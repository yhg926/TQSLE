#ifndef IQSLE_H
#define IQSLE_H

#include <argp.h>


typedef struct indexProperty{
    int fragl;
    int readl;
    int K;
    double lamda;
    char * outpath ;
	char * tref_fa ;
} indexProperty_t;

typedef struct quantProperty{
	int num_remaining_args ;
	double fix_neg;
	int is_fa; // read is in fastq format: 0, in fasta: 1
	char **remaining_args;
	char *indexpath ;
    char *outpath ;
} quantProperty_t;

typedef struct args
{
    struct arg_global* global;
    char* name;
} arg_t;



int index_cmd (indexProperty_t* );
int quant_cmd (quantProperty_t*);
int index_cmd_wrapper (struct argp_state *) ;
int index_cmd_wrapper (struct argp_state *) ;

#endif

