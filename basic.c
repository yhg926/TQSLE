#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "basic.h"
const int Basemap[128] =
{
  [0 ... 127] = DEFAULT,
  ['a'] = 0, ['A'] = 0,
  ['c'] = 1, ['C'] = 1,
  ['g'] = 2, ['G'] = 2,
  ['t'] = 3, ['T'] = 3,
};
const char Mapbase[]={'A','C','G','T'};

int nextPrime(int n){
    int j;
    int tag = 0;

    while (1){
        for(j=2;j<=(int)sqrt(n);j++){
            if(n%j == 0){
                tag = 1;
                break;
            }
        }
        if(tag == 1){
            if(n == 0x7FFFFFFF){
                printf("[ERROR] n exceed 0x7FFFFFFF, Can't find a valid prime\n");
                exit(1);
            }
            n++;
            tag = 0;
        }else{
            return n;
        }
    }
}

void print_bkmer(llong bkmer, int K){
	if (K > 32){printf("K %d is too large, should smaller than 32", K);exit(1); };
	for (int p = 1 ; p <= K ; p++){
		llong temp	= bkmer >> ( ( K - p ) * 2 ) ;
		printf("%c",Mapbase[temp % 4]);	
	}
//	printf("\n");

}


