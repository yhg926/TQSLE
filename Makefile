CC = gcc

CFLAGS = -g -Wall -O2 \
-lrt \
-gdwarf-2 \
-DNPARTITION

SUITESPARSE = /home/yihuiguang/tools/SuiteSparse-5.10.1

LPATH = $(SUITESPARSE)/lib
LIBS = $(LPATH)/libcholmod.so \
$(LPATH)/libamd.so \
$(LPATH)/libcolamd.so \
$(LPATH)/libcamd.so \
$(LPATH)/libccolamd.so \
$(LPATH)/libsuitesparseconfig.so \

CHOLMOD = $(SUITESPARSE)/CHOLMOD
HEADER_DIR = $(CHOLMOD)/Include
CONFIG_HEADER_DIR = $(SUITESPARSE)/SuiteSparse_config

OBJ_DIR = .

BIN_DIR = .

INCLUDES = -I$(HEADER_DIR) \
-I$(CONFIG_HEADER_DIR)

SRCS = $(shell ls *.c)

OBJS = $(SRCS:.c=.o)

OBJS_BUILD = $(shell ls $(OBJ_DIR)/*.o)

APP = iqsle

RM = rm -f

all: $(APP)

$(APP): $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/$(APP) $(OBJS_BUILD) $(LIBS) -lm

%.o: %.c $(HEADER_DIR)/*.h $(CONFIG_HEADER_DIR)/*.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(OBJ_DIR)/$@

clean:
	$(RM) $(OBJS_BUILD) $(APP)
