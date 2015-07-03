
INITIAL  = sedov
HYDRO    = euler
OUTPUT   = ascii

EXE      = SNe


## Set location of HDF5 dir, for use in reading Grackle cooling tables
UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = $(HOME)/bin/yt/yt-x86_64
endif
ifeq ($(UNAME),Darwin)
# assumes MacPorts
H55 = /opt/local
endif

############# Grackle cooling########################################
# See grackle documentation/examples to figure out what's going on here
#
# This'll pollute a lot of the namespace (e.g. DEFINES, CFLAGS, INCLUDES, LIBS)

GRACKLE 			= $(HOME)/bin/grackle
GRACKLE_DIR 		= $(GRACKLE)/src/clib
include $(GRACKLE_DIR)/Make.config.settings

MAKE_CONFIG_MACHINE  = $(GRACKLE_DIR)/Make.config.machine
include $(GRACKLE_DIR)/Make.config.machine

MAKE_CONFIG_OVERRIDE = $(GRACKLE_DIR)/Make.config.override
include $(MAKE_CONFIG_OVERRIDE)
CONFIG_USE_MPI = yes

include $(GRACKLE_DIR)/Make.config.assemble

-include $(GRACKLE_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.grackle/Make.mach.$(CONFIG_MACHINE)

GRACKLE_INCLUDE = -I$(MACH_INSTALL_PREFIX)/include
GRACKLE_LIB = -L$(MACH_INSTALL_PREFIX)/lib -lgrackle


DYLD_LIBRARY_PATH 	= $(GRACKLE)/lib  			
LD_LIBRARY_PATH 	= $(GRACKLE)/lib:$(HOME)/bin/yt/yt-x86_64/lib


############# SUFFIX RULES ########################################


CC = mpicc
FLAGS = -Og -Wall -g 

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lm -lhdf5

OBJ = main.o mpisetup.o profiler.o readpar.o domain.o gridsetup.o geometry.o exchange.o misc.o timestep.o onestep.o riemann.o boundary.o plm.o cooling.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o #report.o

HEADERS = structure.h constants.h

default: $(EXE)

%.o: %.c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(OUTPUT).o : Output/$(OUTPUT).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

$(EXE): $(OBJ) $(HEADERS)
	$(CC) $(LIBS) $(GRACKLE_LIB) $(FLAGS) $(LIB) -o $(EXE) $(OBJ)

clean:
	rm -f *.o $(EXE)
