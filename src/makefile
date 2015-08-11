
INITIAL  = sedov
HYDRO    = euler
OUTPUT   = ascii

EXE      = SNe


############# SYSTEM SPECIFIC ########################################
INC = 
LIB = -lm

ifeq ($(shell uname),Linux)
# General linux options
CC = gcc
CXX = g++

C_FLAGS =  -Wall -g
CXX_FLAGS = -Wall -g -std=c++11

INC += 
LIB += -luuid

INC_BOOST = 
LIB_BOOST =

INC_GSL =
LIB_GSL =

INC_GRACKLE = -I$(HOME)/local/grackle/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(HOME)/local/grackle/lib -lgrackle

DYNLIB = -shared
DYNLIB_EXT = .so
endif

ifeq ($(findstring campusrocks2, $(HOSTNAME)),campusrocks2)
# Campus cluster options
# 	Overwrite general linux options
CC = gcc
CXX = g++

C_FLAGS =  -Wall -g
CXX_FLAGS = -Wall -g -std=c++11

INC += -I$(HOME)/local/libuuid/include
LIB += -L$(HOME)/local/libuuid/lib -luuid

INC_BOOST = 
LIB_BOOST =

INC_GSL =
LIB_GSL =

INC_GRACKLE = -I$(HOME)/local/grackle/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(HOME)/local/grackle/lib -lgrackle

DYNLIB = -shared
DYNLIB_EXT = .so
endif

ifeq ($(shell uname),Darwin)
# Mac options
CC = cc
CXX = c++

C_FLAGS =  -Wall -g
CXX_FLAGS = -Wall -g -std=c++11

INC += 
LIB +=

INC_BOOST = 
LIB_BOOST =

INC_GSL =
LIB_GSL =

INC_GRACKLE = -I$(HOME)/local/grackle/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(HOME)/local/grackle/lib -lgrackle

DYNLIB = -dynamiclib
DYNLIB_EXT = .dylib
endif

INC +=  $(INC_BOOST) $(INC_GSL) $(INC_GRACKLE)
LIB +=  $(LIB_BOOST) $(LIB_GSL) $(LIB_GRACKLE)

LIB_NAME = ../lib/$(EXE)$(DYNLIB_EXT)

############# SYSTEM SPECIFIC ########################################


OBJ = main.o profiler.o readpar.o domain.o gridsetup.o \
      geometry.o misc.o timestep.o riemann.o boundary.o plm.o \
      cooling.o \
      $(INITIAL).o $(OUTPUT).o $(HYDRO).o

HEADERS = structure.h constants.h


############# RULES ########################################


default: all

all: $(EXE)
all: $(LIB_NAME)
all: move_headers

%.o: %.c $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).c $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(OUTPUT).o : Output/$(OUTPUT).c $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Output/$(OUTPUT).c

$(EXE): $(OBJ) $(HEADERS) 
	$(CXX) $(CXX_FLAGS) -o $(EXE) $(OBJ) $(LIB)
	cp $(EXE) ../bin/$(EXE)

$(LIB_NAME): $(OBJ) $(HEADERS)
	$(CXX) $(CXX_FLAGS) -o $(LIB_NAME) $(OBJ) $(LIB)

move_headers: 
	cp $(HEADERS) ../include/

clean:
	rm -f *.o $(EXE)
