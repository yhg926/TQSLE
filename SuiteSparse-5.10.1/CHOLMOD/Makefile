#-------------------------------------------------------------------------------
# CHOLMOD Makefile
#-------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# Compile the C-callable libraries and the Demo programs.
all:
	( cd Lib ; $(MAKE) --jobs=$(JOBS) )

# Compile the C-callable libraries only.
library:
	( cd Lib ; $(MAKE) --jobs=$(JOBS) )

# compile the static libraries only
static:
	( cd Lib    ; $(MAKE) static --jobs=$(JOBS) )

# Remove all files not in the original distribution
purge:
	( cd Lib ; $(MAKE) purge )

# Remove all files not in the original distribution, except keep the 
# compiled libraries.
clean:
	( cd Lib ; $(MAKE) clean )

distclean: purge

ccode: all

# install CHOLMOD
install:
	( cd Lib ; $(MAKE) install )

# uninstall CHOLMOD
uninstall:
	( cd Lib ; $(MAKE) uninstall )
