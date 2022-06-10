#ifndef BASIC
#define BASIC
#include  <limits.h> 

#define _64MASK 0xffffffffffffffffLLU
#define BIT1MASK 0x0000000000000001LLU

#define LMAX 4096
#define LD_FCTR 0.6
#define DEFAULT (-1)
#define HASH_DFLT 0x7FFFFFFF //INT_MAX


typedef unsigned long long int llong;

extern const int Basemap[128];
extern const char Mapbase[];
int nextPrime(int n);

#define H1(K,HASH_SZ) ((K)%(HASH_SZ))
#define H2(K,HASH_SZ) ( 1 + (K) % ( (HASH_SZ) - 1 ) )
#define HASH(K,I,HASH_SZ) ( ( H1(K,HASH_SZ) + I * H2(K,HASH_SZ) ) % HASH_SZ )

void print_bkmer(llong bkmer, int K);

#define SWAP2  0x3333333333333333ULL
#define SWAP4  0x0F0F0F0F0F0F0F0FULL
#define SWAP8  0x00FF00FF00FF00FFULL
#define SWAP16 0x0000FFFF0000FFFFULL
#define SWAP32 0x00000000FFFFFFFFULL

static inline llong crvs64bits(llong n) {

  n = ((n >> 2 ) & SWAP2 ) | ((n & SWAP2 ) << 2 );
  n = ((n >> 4 ) & SWAP4 ) | ((n & SWAP4 ) << 4 );
  n = ((n >> 8 ) & SWAP8 ) | ((n & SWAP8 ) << 8 );
  n = ((n >> 16) & SWAP16) | ((n & SWAP16) << 16);
  n = ((n >> 32) & SWAP32) | ((n & SWAP32) << 32);
  return ~n;
}

typedef struct tref_cat_binary {
    int n; //length of catted ref ( # of nt )
    llong *cat_tref;
    llong *rc_cat_tref; //reverse complete
    int noN; //numer of invalid bases;
    int * Npos; //invalid base positions
    int noT; // numer of tref;
    int *tref_nxt_start; // i+1 trnscrpt start pos, can fast caculate i_start_pos

} tref_cat_binary_t ;


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
            temp +=  ((btref->rc_cat_tref)[nllong - 2 - ind] << (64-(rmd-lshft)) ) >> lshft;
    }

    return (temp);

}

#endif
