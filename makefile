
INITIAL  = sedov
HYDRO    = euler
OUTPUT   = ascii

EXE      = SNe

############# Grackle cooling ########################################
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

############# DEFINITIONS ########################################

ifeq ($(shell uname),Linux)
# General linux options
SYSTEM_INC = 
SYSTEM_LIB = -luuid
ifeq ($(findstring campusrocks2, $(HOSTNAME)),campusrocks2)
# Campus cluster options
# 	Overwrite general linux options
SYSTEM_INC = -I$(HOME)/bin/libuuid/include
SYSTEM_LIB = -L$(HOME)/bin/libuuid/lib -luuid
endif
endif

ifeq ($(shell uname),Darwin)
# Mac options
SYSTEM_INC = 
SYSTEM_LIB = 
endif

CC = mpicc
FLAGS =  -Wall -g 

INC = $(SYSTEM_INC)
LIB = -lm $(SYSTEM_LIB)

OBJ = main.o mpisetup.o profiler.o readpar.o domain.o gridsetup.o geometry.o misc.o timestep.o riemann.o boundary.o plm.o cooling.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o #report.o

HEADERS = structure.h constants.h

############# RULES ########################################


default: $(EXE)

%.o: %.c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INC) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INC) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INC) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INC) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) -c Hydro/$(HYDRO).c

$(OUTPUT).o : Output/$(OUTPUT).c $(HEADERS)
	$(CC) $(DEFINES) $(CFLAGS) $(INC) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS) -c Output/$(OUTPUT).c

$(EXE): $(OBJ) $(HEADERS)
	$(CC) $(FLAGS) -o $(EXE) $(OBJ) $(LIB) $(LIBS) $(GRACKLE_LIB)

clean:
	rm -f *.o $(EXE)
