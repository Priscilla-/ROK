#
#   This is the ROK Code MAKEFILE
#   
#   Last Modification: 26.04.2013 
#

# Fundamental Variables

VERSION = 1.7
NAME = ROK

# Variables Independent of the System
SYSTEM = mac
BASE = $(shell /bin/pwd)
SRCD = $(BASE)/src
EXED = $(BASE)/exe
MAND = $(BASE)/doc/man
INCD = $(BASE)/include


CC = /usr/bin/gcc
CPP = /usr/bin/g++
FC = /usr/bin/g77
INCLUDE = -I$(BASE)/include -I/usr/local/include -I/usr/lib
CFLAGS = -DGNU_MATH $(INCLUDE) -c -W -Wall -O3
CPPFLAGS = -DGNU_MATH $(INCLUDE) -c
#CFLAGS = -I $(INCLUDE) -c -ggdb3 -W -Wall -O3
    ifeq ($(PROFILER),yes)
        LIBRARIES = -pg
        LINK_LIBRARIES = -pg -L/usr/lib -L/usr/local/lib -L/usr/local/ -L/opt/local/lib -L/opt/ioa/software/fftw/3.2.2_double/lib -L/opt/ioa/software -L/opt/ioa/software/gsl/1.15/lib -lgsl -lgslcblas -lm
    else
        LIBRARIES =
        LINK_LIBRARIES = -L/usr/lib -L/usr/local/lib -L/usr/local/ -L/opt/local/lib -L/opt/ioa/software/fftw/3.2.2_double/lib -L/opt/ioa/software -L/opt/ioa/software/gsl/1.15/lib -lgsl -lgslcblas -lm
    endif
#LIBRARIES_DIR = -L/usr/lib -L/usr/local/lib -L/usr/local

# REMOVE COMMAND
REMOVE = /bin/rm


###############################################################################
#
#  VARIABLES TO BE EXPORTED
#
###############################################################################
	
export NAME VERSION BASE SYSTEM EXED CC CPP FC CFLAGS CPPFLAGS LIBRARIES LIBRARIES_DIR LINK_LIBRARIES RM
		

###############################################################################
#
#  TARGETS ASSOCIATED WITH MAINTENANCE OF THE CODE
#
###############################################################################	

.PHONY : clean cleand ccore

info:
	@echo ------------------------------------------------------------------------------------------
	@echo
	@echo   This is the Makefile of the code $(NAME) $(VERSION)
	@echo
	@echo ------------------------------------------------------------------------------------------

clean:
	$(RM) *~
	cd $(INCD); $(RM) *~
	$(MAKE) -C $(SRCD)/main clean

cleand:
	cd data; $(RM) *dat*
	cd data/GnuPlot/; $(RM) -r *

ccore:
	$(RM) core.*



###############################################################################
#
#  ALL CODE TARGETS
#
###############################################################################

all: rok thatsall


###############################################################################
#
#  FROM NOW ON WE SPECIFY TARGETS DIFFERENT FROM THE MAIN EXECUTABLE
#
###############################################################################

rok:
	@echo -----------------------------------------------------------------
	@echo
	@echo   Compiling the $(NAME) code.   Version: $(VERSION)
	@echo
	@echo -----------------------------------------------------------------
	$(MAKE) -C $(SRCD)/main 
			
thatsall:
	@echo
	@echo
	@echo
	@echo
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo
	@echo \ \ [*]  ALL THE $(NAME) EXECUTABLES HAVE BEEN COMPILED!
	@echo
	@echo \ \ \ \ \ \ $(EXED)
	@echo
	@echo \ \ [*]  To run $(NAME): ./exe/rok XXXX
	@echo
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

