
INITIAL  = cluster_SNe
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
LIB += -lstdc++

INC_UUID = 
LIB_UUID =  -luuid

SLUG_DIR = $(HOME)/local/slug2
INC_SLUG = -I$(SLUG_DIR)/src
LIB_SLUG = -L$(SLUG_DIR)/src -lslug

INC_BOOST =  -I$(HOME)/local/boost_1_59_0/include
LIB_BOOST =  -L$(HOME)/local/boost_1_59_0/lib -lboost_system -lboost_filesystem

INC_GSL = -I$(HOME)/local/gsl-1.16/include
LIB_GSL = -L$(HOME)/local/gsl-1.16/lib

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR = $(HOME)/local/grackle
INC_GRACKLE = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(GRACKLE_DIR)/lib -lgrackle

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

INC += 
LIB += -lstdc++

INC_UUID = -I$(HOME)/local/libuuid/include
LIB_UUID = -L$(HOME)/local/libuuid/lib -luuid

SLUG_DIR = $(HOME)/local/slug2
INC_SLUG = -I$(SLUG_DIR)/src
LIB_SLUG = -L$(SLUG_DIR)/src -lslug

INC_BOOST =  -I$(HOME)/local/boost_1_59_0/include
LIB_BOOST =  -L$(HOME)/local/boost_1_59_0/lib -lboost_system -lboost_filesystem

INC_GSL = -I$(HOME)/local/gsl-1.16/include
LIB_GSL = -L$(HOME)/local/gsl-1.16/lib

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR = $(HOME)/local/grackle
INC_GRACKLE = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(GRACKLE_DIR)/lib -lgrackle

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

INC_UUID = 
LIB_UUID = 


SLUG_DIR = $(HOME)/local/slug2
INC_SLUG = -I$(SLUG_DIR)/src
LIB_SLUG = -L$(SLUG_DIR)/src -lslug

INC_BOOST = 
LIB_BOOST = -lboost_system -lboost_filesystem

INC_GSL =
LIB_GSL =

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR = $(HOME)/local/grackle
INC_GRACKLE = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(GRACKLE_DIR)/lib -lgrackle

DYNLIB = -dynamiclib
DYNLIB_EXT = .dylib
endif

INC +=  $(INC_UUID) $(INC_SLUG) $(INC_BOOST) $(INC_GSL) $(INC_GRACKLE)
LIB +=  $(LIB_UUID) $(LIB_SLUG) $(LIB_BOOST) $(LIB_GSL) $(LIB_GRACKLE)


INC_DIR = ../include/
LIB_DIR = ../lib/
LIB_NAME = $(LIB_DIR)lib$(EXE)$(DYNLIB_EXT)

DEFINES += -DGRACKLE_DIR=\"$(GRACKLE_DIR)\"
DEFINES += -DSLUG_DIR=\"$(SLUG_DIR)\"

############# SYSTEM SPECIFIC ########################################


OBJ = main.o profiler.o readpar.o domain.o gridsetup.o \
      geometry.o misc.o timestep.o riemann.o boundary.o plm.o \
      cooling.o blast.o \
      $(INITIAL).o $(OUTPUT).o $(HYDRO).o

HEADERS = structure.H constants.H blast.H domain.H readpar.H

############# RULES ########################################


default: all

all: $(LIB_NAME)
all: $(EXE)
all: move_headers

%.o: %.c $(HEADERS)
	$(CC) $(DEFINES) $(C_FLAGS) $(INC) -c $<

%.o: %.cxx $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c $<

$(TIMESTEP).o: Timestep/$(TIMESTEP).cxx $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Timestep/$(TIMESTEP).cxx

$(INITIAL).o : Initial/$(INITIAL).cxx $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Initial/$(INITIAL).cxx

$(HYDRO).o : Hydro/$(HYDRO).cxx $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Hydro/$(HYDRO).cxx

$(OUTPUT).o : Output/$(OUTPUT).cxx $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c Output/$(OUTPUT).cxx

$(EXE): $(OBJ) $(HEADERS) $(LIB_NAME)
	$(CXX) $(CXX_FLAGS) -o $(EXE) -L../lib -lSNe $(LIB)
	cp $(EXE) ../bin/$(EXE)

$(LIB_NAME): $(OBJ) $(HEADERS)
	mkdir -p $(LIB_DIR)
	$(CXX) $(CXX_FLAGS) -o $(LIB_NAME) $(OBJ) $(LIB) $(DYNLIB)

move_headers: 
	mkdir -p $(INC_DIR)
	cp $(HEADERS) $(INC_DIR)

clean:
	rm -f *.o $(EXE) $(LIB_NAME) ../bin/$(EXE)
