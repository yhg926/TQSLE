CC = gcc

CFLAGS =-std=gnu99 -g -Wall -O2 \
-lrt \
-gdwarf-2

OPENBLAS = $(CURDIR)/OpenBLAS-0.3.18
SUITESPARSE = $(CURDIR)/SuiteSparse-5.10.1

LPATH = $(SUITESPARSE)/lib
LIBS = -L$(LPATH) -Wl,-rpath=$(LPATH) -lcholmod -lamd -lcolamd -lcamd -lccolamd -lsuitesparseconfig 

CHOLMOD = $(SUITESPARSE)/CHOLMOD
CHOLMOD_HEADER_DIR = $(CHOLMOD)/Include
HEADER_DIR = $(SUITESPARSE)/include #$(CHOLMOD)/Include
CONFIG_HEADER_DIR = $(SUITESPARSE)/SuiteSparse_config

OBJ_DIR = .

BIN_DIR = .

INCLUDES = -I$(HEADER_DIR) \
-I$(CHOLMOD_HEADER_DIR) \
-I$(CONFIG_HEADER_DIR)

SRCS = $(shell ls *.c)

OBJS = $(SRCS:.c=.o)

OBJS_BUILD = $(shell ls $(OBJ_DIR)/*.o)

APP = iqsle

RM = rm -f

CP = cp 

all: need $(APP) install clean

need: 
	( cd $(OPENBLAS) && make && make PREFIX=$(OPENBLAS) install )
	( cd $(SUITESPARSE) && make BLAS="-L$(OPENBLAS)/lib -Wl,-rpath=$(OPENBLAS)/lib -lopenblas" )

$(APP): $(OBJS) 
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$(APP) $(OBJS_BUILD) $(LIBS) -lm

%.o: %.c 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(OBJ_DIR)/$@

clean:
	$(RM) $(OBJS_BUILD) $(APP)

install:
	mkdir -p $(CURDIR)/bin/
	$(CP) $(APP) $(CURDIR)/bin/

uninstall:
	( cd $(OPENBLAS) && make clean)
	( cd $(SUITESPARSE) && make uninstall)
	$(RM) $(CURDIR)/bin/*
