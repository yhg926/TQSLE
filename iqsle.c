#include <argp.h>
#include <string.h>
#include <stdlib.h>
#include <error.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "basic.h"
#include "iqsle.h"
#include "index_core.h"


const char version[] = "IQSLE v1.0";

//index command args-------- 

static struct argp_option indexOpt[] = {
    {"kmerLen", 'k', "INT", 0, "Length of kmer [31].\v"},
	{"readLen", 'r', "INT", 0, "Length of read [100]\v"},
	{"fragLen", 'f', "INT", 0, "Length of fragment [300]\v"},
    {"lamda", 'p', "DOUBLE", 0, "for solving (A'A + lamda*I)x = A'b [10].\v"},
    {"output", 'o', "STRING", 0, "Index output path.\v"},
	{"input", 'i', "STRING", 0, "transcript reference (.fasta).\v"},
    { 0 }
};

static char indexArgsDoc[] = "-i reference.fasta";

static indexProperty_t indexProperty = {
	300, //frag len 
	100, //read len
	 31, //kmer K
	 10, //lamdai
	 "./", //output path
	NULL
	};

const char indexDoc[] = "\nWelcome to use IQSLE v1.0\n";

static error_t indexParse(int key, char* arg, struct argp_state* state) {

    switch(key)
    {
        case 'k':
            {
                indexProperty.K = atoi(arg);
                break;
            }
        case 'r':
            {
                indexProperty.readl = atoi(arg);
                break;
            }
        case 'f':
            {
                indexProperty.fragl = atoi(arg);
                break;
            }
        case 'p':
            {
                indexProperty.lamda = atof(arg);
                break;
            }
        case 'o':
            {
                indexProperty.outpath = arg ;

				struct stat st = {0};
				if (stat(indexProperty.outpath, &st) == -1) {
    				mkdir(indexProperty.outpath, 0700);
				}
                break;
            }
		case 'i':
			{
				indexProperty.tref_fa = arg ;
				break;
			}

        case ARGP_KEY_NO_ARGS:
            {
				if (  indexProperty.tref_fa  != NULL )  break;
                printf("\v");
                argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
				printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_LONG);	
				printf("\v");
                exit (0);
							
            }
            break;

        default:
            return ARGP_ERR_UNKNOWN;

    }
    return 0;
}

static struct argp indexArgp = {indexOpt, indexParse, indexArgsDoc, indexDoc};

int index_cmd_wrapper (struct argp_state *state){

    arg_t index_info = { .global  = state->input, };   
    int argc = state->argc - 1;
    char **argv = (state->argv)+1;

    argv[0] = malloc(strlen(state->argv[0]) + strlen(" index") + 1);
    sprintf(argv[0], "%s index", state->argv[0]);

    argp_parse (&indexArgp, argc, argv, ARGP_IN_ORDER, &argc, &index_info);

    return index_cmd( &indexProperty);
}

//-------------quant command args---


static struct argp_option quantOpt[] = {
	{"indexPath", 'i', "STRING", 0, "index path.\v"},
    {"output", 'o', "STRING", 0, "output path.\v"},
    { 0 }
};

static char quantArgsDoc[] = "*.fastq [...]";

static quantProperty_t quantProperty = { 
	0,//int num_remaining_args ;
	NULL, //char **remaining_args;	
	"./" ,//char *indexpath ; 
 	"./"  //char *outpath ;
};

static char quantDoc[] = "\nWelcome to use rabbit v0.1\n";

static error_t quantParse(int key, char* arg, struct argp_state* state) {

    switch(key)
    {
        case 'o':
            {
                quantProperty.outpath = arg ;

                struct stat st = {0};
                if (stat(quantProperty.outpath, &st) == -1) {
                    mkdir(quantProperty.outpath, 0700);
                }
                break;
            }
        case 'i':
            {
                quantProperty.indexpath = arg;
                break;
            }
	  case ARGP_KEY_ARGS:
		{
    		quantProperty.num_remaining_args = state->argc - state->next;
        	quantProperty.remaining_args = state->argv + state->next;
   			break;
		}
        case ARGP_KEY_NO_ARGS:
            {
                printf("\v");
                argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
                printf("\v");
                argp_state_help(state,stdout,ARGP_HELP_LONG);
                printf("\v");				
                exit(0);
            }

        default:
            return ARGP_ERR_UNKNOWN;
	}
    return 0;
}

static struct argp quantArgp = {quantOpt, quantParse, quantArgsDoc, quantDoc};

int quant_cmd_wrapper(struct argp_state *state){	

	arg_t quant_info = { .global  = state->input, };
    int argc = state->argc - 1;
    char **argv = (state->argv)+1;

    argv[0] = malloc(strlen(state->argv[0]) + strlen(" quant") + 1);
    sprintf(argv[0], "%s quant", state->argv[0]);

    argp_parse (&quantArgp, argc, argv, ARGP_IN_ORDER, &argc, &quant_info);

    return  quant_cmd (&quantProperty);
}



//-------------main--------

static struct argp_option mainOpt[] = {
    {"version", 'v', 0, OPTION_NO_USAGE, 0},
    { 0 }
};

static char mainArgsDoc[] = "ARG [STRING...]";

static char mainDoc[] =
"\n"
    "  Command:\n"
"\n"
    "    index\t\tmake transcripts index\n"
"\n"
    "    quant\t\tisoforms expression level quantification\n"
"\n";

static error_t mainParse(int key, char* arg, struct argp_state* state) {

    if(key ==  'v'){
            printf("\v\t%s\n",version);
            return (0);
    }else if (key == ARGP_KEY_ARG){
            if (strcmp(arg, "index") == 0){
                index_cmd_wrapper(state);
                return(1);
            }else if(strcmp(arg, "quant") == 0){
                quant_cmd_wrapper(state);
                return(2);
            }else{
                argp_error(state, "%s is not a valid command", arg);
            }
    }
	else if (key == ARGP_KEY_NO_ARGS){
            printf("\v");
            argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
            printf("%s\n", mainDoc);
            printf("\v");
			return -1;
    }
    return 0;
}

static struct argp argp = {mainOpt, mainParse, mainArgsDoc, mainDoc, 0, 0, 0};

int main(int argc, char **argv){

    argp_parse (&argp, argc, argv, ARGP_IN_ORDER, &argc, 0);

    return 0;
}


