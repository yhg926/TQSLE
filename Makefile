CC = gcc

OFLAGS =-std=gnu99 -g -Wall -O2 -fopenmp
CFLAGS =-std=gnu99 -g -Wall -O2 \
-lrt \
-gdwarf-2

OMP_CMD =-fopenmp


BLAS_OPT =  

MKL_LIB =
MKL_IOMP5 =
BLAS_LIB =

OPENBLAS = $(CURDIR)/OpenBLAS-0.3.18
SUITESPARSE = $(CURDIR)/SuiteSparse-5.10.1

ifdef MKL_LIB
	ifdef MKL_IOMP5
		BLAS_OPT = mkl
		OMP_CMD =  -L$(MKL_IOMP5) -Wl,-rpath=$(MKL_IOMP5) -liomp5
		BLAS_CMD = -L$(MKL_LIB) -Wl,-rpath=$(MKL_LIB) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread $(OMP_CMD) -lpthread -lm
	endif
endif


ifndef BLAS_OPT
	ifdef BLAS_LIB
		BLAS_OPT = blas
		BLAS_CMD = -L$(BLAS_LIB) -Wl,-rpath=$(BLAS_LIB) -lopenblas
	endif
endif

ifndef BLAS_OPT
	BLAS_OPT = local
	BLAS_LIB = $(OPENBLAS)/lib
	BLAS_CMD = -L$(BLAS_LIB) -Wl,-rpath=$(BLAS_LIB) -lopenblas
endif


LPATH = $(SUITESPARSE)/lib
LIBS = -L$(LPATH) -Wl,-rpath=$(LPATH) -lcholmod -lamd -lcolamd -lcamd -lccolamd -lsuitesparseconfig 

CHOLMOD = $(SUITESPARSE)/CHOLMOD
CHOLMOD_HEADER_DIR = $(CHOLMOD)/Include
HEADER_DIR = $(SUITESPARSE)/include #$(CHOLMOD)/Include
CONFIG_HEADER_DIR = $(SUITESPARSE)/SuiteSparse_config

OBJ_DIR = $(CURDIR)/bin

BIN_DIR = $(CURDIR)/bin

INCLUDES = -I$(HEADER_DIR) \
-I$(CHOLMOD_HEADER_DIR) \
-I$(CONFIG_HEADER_DIR)

SRCS = $(shell ls *.c)

OBJS = $(SRCS:.c=.o)

OBJS_BUILD = $(shell ls $(OBJ_DIR)/*.o)

APP = iqsle

RM = rm -f

CP = cp 

MKDIR = mkdir -p

all: need $(APP) clean

need:  
ifeq ($(BLAS_OPT),local)
	( cd $(OPENBLAS) && make && make PREFIX=$(OPENBLAS) install )
endif
	( cd $(SUITESPARSE) && make BLAS="$(BLAS_CMD)" )
	$(MKDIR) $(OBJ_DIR) $(BIN_DIR) 

IOMP5_DIR =
IOMP5 = $(shell ldd $(LPATH)/libcholmod.so | grep 'iomp5' | cut -f3 -d" " )

ifneq ($(IOMP5),) #use ifdef cause error
	IOMP5_DIR = $(shell dirname $(IOMP5))
endif

ifneq ($(IOMP5_DIR),)
	OMP_CMD = -L$(IOMP5_DIR) -Wl,-rpath=$(IOMP5_DIR) -liomp5
endif


$(APP): $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$(APP) $(OBJS_BUILD) $(LIBS) $(OMP_CMD) -lm

%.o: %.c 
	$(CC) $(OFLAGS) $(INCLUDES) -c $< -o $(OBJ_DIR)/$@ 

clean:
	$(RM) $(OBJS_BUILD)

uninstall:
	( cd $(SUITESPARSE) && make uninstall)
	$(RM) $(BIN_DIR)/*

check:
	( ldd $(BIN_DIR)/$(APP) | grep 'mkl\|blas\|omp\|not found')
