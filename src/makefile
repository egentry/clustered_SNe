
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

SLUG_DIR  = $(HOME)/local/slug2
INC_SLUG  = -I$(SLUG_DIR)/src
LIB_SLUG  = -L$(SLUG_DIR)/src -lslug
LIB_SLUG += -Wl,-rpath,$(SLUG_DIR)/src

BOOST_DIR  = $(HOME)/local/boost_1_59_0
INC_BOOST  = -I$(BOOST_DIR)/include
LIB_BOOST  = -L$(BOOST_DIR)/lib -lboost_system -lboost_filesystem
LIB_BOOST += -Wl,-rpath,$(BOOST_DIR)/lib

GSL_DIR  = $(HOME)/local/gsl-1.16
INC_GSL  = -I$(GSL_DIR)/include
LIB_GSL  = -L$(GSL_DIR)/lib 
LIB_GSL += -Wl,-rpath,$(HOME)/local/gsl-1.16/lib

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR  =  $(HOME)/local/grackle
INC_GRACKLE  = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE  = -L$(GRACKLE_DIR)/lib -lgrackle
LIB_GRACKLE += -Wl,-rpath,$(GRACKLE_DIR)/lib

DYNLIB = -shared
DYNLIB_EXT = .so
endif

ifeq ($(findstring campusrocks2, $(HOSTNAME)),campusrocks2)
# Campus cluster options
# 	Overwrite general linux options
CC = gcc
CXX = g++

C_FLAGS =  -Wall -g
CXX_FLAGS = -Wall -g -std=c++11 -fPIC

INC += 
LIB += -lstdc++

UUID_DIR = $(HOME)/local/libuuid
INC_UUID = -I$(UUID_DIR)/include
LIB_UUID = -L$(UUID_DIR)/lib -luuid

SLUG_DIR  = $(HOME)/local/slug2
INC_SLUG  = -I$(SLUG_DIR)/src
LIB_SLUG  = -L$(SLUG_DIR)/src -lslug
LIB_SLUG += -Wl,-rpath,$(SLUG_DIR)/src

BOOST_DIR  = $(HOME)/local/boost_1_59_0
INC_BOOST  = -I$(BOOST_DIR)/include
LIB_BOOST  = -L$(BOOST_DIR)/lib -lboost_system -lboost_filesystem
LIB_BOOST += -Wl,-rpath,$(BOOST_DIR)/lib

GSL_DIR  = $(HOME)/local/gsl-1.16
INC_GSL  = -I$(GSL_DIR)/include
LIB_GSL  = -L$(GSL_DIR)/lib 
LIB_GSL += -Wl,-rpath,$(HOME)/local/gsl-1.16/lib

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR  =  $(HOME)/local/grackle
INC_GRACKLE  = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE  = -L$(GRACKLE_DIR)/lib -lgrackle
LIB_GRACKLE += -Wl,-rpath,$(GRACKLE_DIR)/lib

#cfitsio an optional dependency of slug2
ETC_RPATHS  = -Wl,-rpath,$(HOME)/local/cfitsio/lib
ETC_RPATHS += -Wl,-rpath,$(HOME)/local/gcc-5.2.0/lib64
LIB += $(ETC_RPATHS)

DYNLIB = -shared
DYNLIB_EXT = .so
endif

ifeq ($(shell uname),Darwin)
# Mac options
CC = cc
CXX = c++

C_FLAGS =  -Wall -g
CXX_FLAGS = -Wall -g -std=c++11 -fPIC


INC += 
LIB += 

INC_UUID = 
LIB_UUID = 


SLUG_DIR = $(HOME)/local/slug2
INC_SLUG = -I$(SLUG_DIR)/src
LIB_SLUG = -L$(SLUG_DIR)/src -lslug -Wl,-rpath,$(SLUG_DIR)/src

INC_BOOST = 
LIB_BOOST = -lboost_system -lboost_filesystem

INC_GSL =
LIB_GSL =

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR = $(HOME)/local/grackle
INC_GRACKLE = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(GRACKLE_DIR)/lib -lgrackle -Wl,-rpath,$(GRACKLE_DIR)/lib

DYNLIB = -dynamiclib
DYNLIB_EXT = .dylib
LIB_FLAGS = -install_name @rpath/$(EXE)$(DYNLIB_EXT)
endif

INC +=  $(INC_UUID) $(INC_SLUG) $(INC_BOOST) $(INC_GSL) $(INC_GRACKLE)
LIB +=  $(LIB_UUID) $(LIB_SLUG) $(LIB_BOOST) $(LIB_GSL) $(LIB_GRACKLE)

PREFIX = ..

INC_DIR = $(PREFIX)/include
LIB_DIR = $(PREFIX)/lib
LIB_NAME = $(LIB_DIR)/$(EXE)$(DYNLIB_EXT)

DEFINES += -DGRACKLE_DIR=\"$(GRACKLE_DIR)\"
DEFINES += -DSLUG_DIR=\"$(SLUG_DIR)\"

############# SYSTEM SPECIFIC ########################################


OBJ = main.o profiler.o readpar.o domain.o \
      geometry.o misc.o timestep.o riemann.o boundary.o plm.o \
      cooling.o blast.o \
      Initial/initial_conditions.o \
      Initial/restart_ICs.o Initial/chevalier_ICs.o \
      Initial/cluster_SNe_ICs.o Initial/ejecta_ICs.o \
      Initial/isentropic_ICs.o Initial/messy_ICs.o \
      Initial/shocktube_ICs.o Initial/Thornton_parameter_study_ICs.o \
      Initial/uniform_ICs.o \
      Output/ascii.o Hydro/euler.o

HEADERS = 	blast.H boundary.H constants.H cooling.H domain.H \
			geometry.H misc.H plm.H profiler.H readpar.H \
			structure.H timestep.H \
	      	Output/ascii.H Hydro/euler.H \
			Initial/Thornton_parameter_study_ICs.H \
			Initial/chevalier_ICs.H Initial/cluster_SNe_ICs.H \
			Initial/ejecta_ICs.H Initial/initial_conditions.H \
			Initial/isentropic_ICs.H Initial/messy_ICs.H \
			Initial/restart_ICs.H Initial/shocktube_ICs.H \
			Initial/uniform_ICs.H


############# RULES ########################################

default: all 

all: $(LIB_NAME)
all: $(EXE)
all: move_headers


%.o: %.cxx $(HEADERS)
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -c $< -o $@

$(EXE): $(OBJ) $(HEADERS) $(LIB_NAME)
	$(CXX) $(CXX_FLAGS) -o $(EXE) $(OBJ) $(LIB)
	cp $(EXE) ../bin/$(EXE)

exe_from_lib: $(OBJ) $(HEADERS) $(LIB_NAME)
	$(CXX) $(CXX_FLAGS) -o $(EXE) -L../lib -lSNe -Wl,-rpath,$(PWD)/../lib $(LIB)
	cp $(EXE) ../bin/$(EXE)

$(LIB_NAME): $(OBJ) $(HEADERS)
	mkdir -p $(LIB_DIR)
	$(CXX) $(CXX_FLAGS) -o $(LIB_NAME) $(OBJ) $(LIB) $(DYNLIB) $(LIB_FLAGS) -LInitial/

move_headers: 
	@ mkdir -p $(INC_DIR)
	@ mkdir -p $(INC_DIR)/Initial
	@ mkdir -p $(INC_DIR)/Hydro
	@ mkdir -p $(INC_DIR)/Output
	@ mkdir -p $(INC_DIR)/Initial
	@ echo "Moving headers into " $(INC_DIR)
	@ for HEADER in $(HEADERS); \
	do \
		cp -f $$HEADER $(INC_DIR)/$$HEADER ; \
	done

clean:
	rm -f *.o $(EXE) $(LIB_NAME) ../bin/$(EXE)
	rm -f Hydro/*.o Output/*.o
