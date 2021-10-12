#-------------------------------------------------------------------------------
# Makefile for all SuiteSparse packages
#-------------------------------------------------------------------------------

SUITESPARSE = $(CURDIR)
export SUITESPARSE

default: go

include SuiteSparse_config/SuiteSparse_config.mk

go:
	( cd SuiteSparse_config && $(MAKE) )
	( cd AMD && $(MAKE) )
	( cd CAMD && $(MAKE) )
	( cd CCOLAMD && $(MAKE) )
	( cd COLAMD && $(MAKE) )
	( cd CHOLMOD && $(MAKE) )

# install all packages in SuiteSparse/lib and SuiteSparse/include.  Use the
# following command to install in /usr/local/lib and /usr/local/include:
#       sudo make install INSTALL=/usr/local
# See SuiteSparse/README.md for more details.
# (note that CSparse is not installed; CXSparse is installed instead)
install:
	( cd SuiteSparse_config && $(MAKE) install )
	( cd AMD && $(MAKE) install )
	( cd CAMD && $(MAKE) install )
	( cd CCOLAMD && $(MAKE) install )
	( cd COLAMD && $(MAKE) install )
	( cd CHOLMOD && $(MAKE) install )

# uninstall all packages
uninstall:
	$(RM) $(INSTALL_DOC)/SuiteSparse_README.md
	( cd SuiteSparse_config && $(MAKE) uninstall )
	( cd AMD && $(MAKE) uninstall )
	( cd CAMD && $(MAKE) uninstall )
	( cd COLAMD && $(MAKE) uninstall )
	( cd CCOLAMD && $(MAKE) uninstall )
	( cd CHOLMOD && $(MAKE) uninstall )
	$(RM) -r $(INSTALL_DOC)

library:
	( cd SuiteSparse_config && $(MAKE) )
	( cd AMD && $(MAKE) library )
	( cd CAMD && $(MAKE) library )
	( cd CCOLAMD && $(MAKE) library )
	( cd COLAMD && $(MAKE) library )
	( cd CHOLMOD && $(MAKE) library )

static:
	( cd SuiteSparse_config && $(MAKE) static )
	( cd AMD && $(MAKE) static )
	( cd CAMD && $(MAKE) static )
	( cd CCOLAMD && $(MAKE) static )
	( cd COLAMD && $(MAKE) static )
	( cd CHOLMOD && $(MAKE) static )

# Remove all files not in the original distribution
purge:
	- ( cd SuiteSparse_config && $(MAKE) purge )
	- ( cd AMD && $(MAKE) purge )
	- ( cd CAMD && $(MAKE) purge )
	- ( cd COLAMD && $(MAKE) purge )
	- ( cd CCOLAMD && $(MAKE) purge )
	- ( cd CHOLMOD && $(MAKE) purge )
	- $(RM) -r include/* bin/* lib/* share/*

# Remove all files not in the original distribution, but keep the libraries
clean:
	- ( cd SuiteSparse_config && $(MAKE) clean )
	- ( cd AMD && $(MAKE) clean )
	- ( cd CAMD && $(MAKE) clean )
	- ( cd COLAMD && $(MAKE) clean )
	- ( cd CCOLAMD && $(MAKE) clean )
	- ( cd CHOLMOD && $(MAKE) clean )

distclean: purge


