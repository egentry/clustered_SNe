
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
LIB += -luuid

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

INC +=  $(INC_SLUG) $(INC_BOOST) $(INC_GSL) $(INC_GRACKLE)
LIB +=  $(LIB_SLUG) $(LIB_BOOST) $(LIB_GSL) $(LIB_GRACKLE)

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

HEADERS = structure.H constants.H blast.H domain.H

############# RULES ########################################


default: all

all: $(EXE)
all: $(LIB_NAME)
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

$(EXE): $(OBJ) $(HEADERS) 
	$(CXX) $(CXX_FLAGS) -o $(EXE) $(OBJ) $(LIB)
	cp $(EXE) ../bin/$(EXE)

$(LIB_NAME): $(OBJ) $(HEADERS) $(EXE)
	mkdir -p $(LIB_DIR)
	$(CXX) $(CXX_FLAGS) -o $(LIB_NAME) $(OBJ) $(LIB)

move_headers: 
	mkdir -p $(INC_DIR)
	cp $(HEADERS) $(INC_DIR)

clean:
	rm -f *.o $(EXE)
