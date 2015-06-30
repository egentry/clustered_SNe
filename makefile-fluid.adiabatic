CC=gcc

FLAGS 		= -lm 
FLAGS_DEBUG = -ggdb -Wall -lm -Og -p

EXE_FILE 	= fluid.adiabatic

OBJS 		= fluid.adiabatic.o \
				cooling.o \
				ICs.o \
				grid.o


############# Grackle cooling########################################
# See grackle documentation/examples to figure out what's going on here
#
# This'll pollute the namespace (e.g. DEFINES, CFLAGS, INCLUDES, LIBS)

GRACKLE 			= $(HOME)/bin/grackle
GRACKLE_DIR 		= $(GRACKLE)/src/clib
include $(GRACKLE_DIR)/Make.config.settings

MAKE_CONFIG_MACHINE  = $(GRACKLE_DIR)/Make.config.machine
include $(GRACKLE_DIR)/Make.config.machine

MAKE_CONFIG_OVERRIDE = $(GRACKLE_DIR)/Make.config.override
include $(MAKE_CONFIG_OVERRIDE)
CONFIG_USE_MPI = no

include $(GRACKLE_DIR)/Make.config.assemble

-include $(GRACKLE_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.grackle/Make.mach.$(CONFIG_MACHINE)

GRACKLE_INCLUDE = -I$(MACH_INSTALL_PREFIX)/include
GRACKLE_LIB = -L$(MACH_INSTALL_PREFIX)/lib -lgrackle


DYLD_LIBRARY_PATH 	= $(GRACKLE)/lib  			
LD_LIBRARY_PATH 	= $(GRACKLE)/lib:$(HOME)/bin/yt/yt-x86_64/lib

############# SUFFIX RULES ########################################

$(EXE_FILE) : $(OBJS)
	@echo "linking..."
	$(CC)  $(OBJS) -o $(EXE_FILE) $(LIBS) $(GRACKLE_LIB) $(FLAGS)

.SUFFIXES :             # clear all defaults
.SUFFIXES : .c .o       # and replace

%.o: %.c %.h constants.h
	$(CC) -c $< $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $(FLAGS)

############# MAKE RULES ########################################

clean:
	@rm -f *.o
	@rm -f $(EXE_FILE)


debug: FLAGS=$(FLAGS_DEBUG)
debug: clean     				# is there a better way to force .o files to be gdb friendly?
debug: $(EXE_FILE)

test:
	@echo DEFINE = $(DEFINE)
	@echo CFLAGS = $(CFLAGS)
	@echo INCLUDES = $(INCLUDES)
	@echo LIBS = $(LIBS)



############# DEPENDENCIES ########################################

fluid.adiabatic.o : cooling.o ICs.o grid.o


