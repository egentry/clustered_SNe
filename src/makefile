
EXE_NAME = SNe

PREFIX = ..

INC_DIR = $(PREFIX)/include
LIB_DIR = $(PREFIX)/lib
BIN_DIR = $(PREFIX)/bin

############# SYSTEM SPECIFIC ########################################
INC = 
LIB = -lm

ifeq ($(shell uname),Linux)
# General linux options
CXX = g++

CXX_FLAGS = -Wall -std=c++11 -fPIC 
CXX_FLAGS_DEBUG = -g -Og
CXX_FLAGS_OPT = -O3 -DNDEBUG

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

ifeq ($(findstring campusrocks, $(HOSTNAME)),campusrocks)
# Campus cluster options
# 	(Overwrite general linux options)
CXX = g++

CXX_FLAGS = -Wall -std=c++11 -fPIC
CXX_FLAGS_DEBUG = -g -Og
CXX_FLAGS_OPT = -O3 -DNDEBUG

INC += 
LIB += -lstdc++

UUID_DIR = $(HOME)/miniconda3
INC_UUID = -I$(UUID_DIR)/include
LIB_UUID = -L$(UUID_DIR)/lib -luuid

SLUG_DIR  = $(HOME)/local/slug2
INC_SLUG  = -I$(SLUG_DIR)/src
LIB_SLUG  = -L$(SLUG_DIR)/src -lslug
LIB_SLUG += -Wl,-rpath,$(SLUG_DIR)/src

BOOST_DIR  = $(HOME)/local/boost_1_60_0
INC_BOOST  = -I$(BOOST_DIR)/include
LIB_BOOST  = -L$(BOOST_DIR)/lib -lboost_system -lboost_filesystem
LIB_BOOST += -Wl,-rpath,$(BOOST_DIR)/lib

GSL_DIR  = $(HOME)/local/gsl-2.1
INC_GSL  = -I$(GSL_DIR)/include
LIB_GSL  = -L$(GSL_DIR)/lib 
LIB_GSL += -Wl,-rpath,$(GSL_DIR)/lib

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR  =  $(HOME)/local/grackle
INC_GRACKLE  = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE  = -L$(GRACKLE_DIR)/lib -lgrackle
LIB_GRACKLE += -Wl,-rpath,$(GRACKLE_DIR)/lib

#cfitsio an optional dependency of slug2
ETC_RPATHS  = -Wl,-rpath,$(HOME)/local/cfitsio/lib
ETC_RPATHS += -Wl,-rpath,$(HOME)/local/gcc-5.3.0/lib64
LIB += $(ETC_RPATHS)

DYNLIB = -shared
DYNLIB_EXT = .so
endif

ifneq ( ,$(filter $(HOSTNAME), hyades.ucsc.edu eudora.ucsc.edu))
# Campus cluster options
# 	(Overwrite general linux options)
CXX = g++

CXX_FLAGS = -Wall -std=c++11 -fPIC
CXX_FLAGS_DEBUG = -g -Og
CXX_FLAGS_OPT = -O3 -DNDEBUG

INC += 
LIB += -lstdc++

UUID_DIR = $(HOME)/local/miniconda3/envs/hdf
INC_UUID = -I$(UUID_DIR)/include 
LIB_UUID = -L$(UUID_DIR)/lib -luuid

SLUG_DIR  = $(HOME)/local/slug2
INC_SLUG  = -I$(SLUG_DIR)/src
LIB_SLUG  = -L$(SLUG_DIR)/src -lslug
LIB_SLUG += -Wl,-rpath,$(SLUG_DIR)/src

BOOST_DIR  = $(HOME)/local/boost_1_60_0
INC_BOOST  = -I$(BOOST_DIR)/include
LIB_BOOST  = -L$(BOOST_DIR)/lib -lboost_system -lboost_filesystem
LIB_BOOST += -Wl,-rpath,$(BOOST_DIR)/lib

GSL_DIR  = $(HOME)/local/gsl-2.1
INC_GSL  = -I$(GSL_DIR)/include
LIB_GSL  = -L$(GSL_DIR)/lib 
LIB_GSL += -Wl,-rpath,$(GSL_DIR)/lib

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR  =  $(HOME)/local/grackle
INC_GRACKLE  = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE  = -L$(GRACKLE_DIR)/lib -lgrackle
LIB_GRACKLE += -Wl,-rpath,$(GRACKLE_DIR)/lib

#cfitsio an optional dependency of slug2
ETC_RPATHS  = -Wl,-rpath,$(HOME)/local/cfitsio/lib
# ETC_RPATHS += -Wl,-rpath,$(HOME)/local/gcc-5.3.0/lib64
LIB += $(ETC_RPATHS)

DYNLIB = -shared
DYNLIB_EXT = .so
endif

ifeq ($(shell uname),Darwin)
# Mac options
CXX = /usr/bin/clang++

CXX_FLAGS = -Wall -std=c++11 -fPIC
CXX_FLAGS_DEBUG = -g -O1
CXX_FLAGS_OPT = -O3 -DNDEBUG 


INC += 
LIB += 

INC_UUID = 
LIB_UUID = 


SLUG_DIR = $(HOME)/local/slug2
INC_SLUG = -I$(SLUG_DIR)/src
LIB_SLUG = -L$(SLUG_DIR)/src -lslug -Wl,-rpath,$(SLUG_DIR)/src

INC_BOOST = -I/usr/local/include
LIB_BOOST = -L/usr/local/lib -lboost_system -lboost_filesystem

INC_GSL =
LIB_GSL =

LAPACKE_DIR = /usr/local/opt/lapack
INC_LAPACKE = -I$(LAPACKE_DIR)/include
LIB_LAPACKE = -L$(LAPACKE_DIR)/lib -llapacke

#GRACKLE_DIR is needed for a define statement within cooling
GRACKLE_DIR = $(HOME)/local/grackle
INC_GRACKLE = -I$(GRACKLE_DIR)/include -DCONFIG_BFLOAT_8
LIB_GRACKLE = -L$(GRACKLE_DIR)/lib -lgrackle -Wl,-rpath,$(GRACKLE_DIR)/lib

DYNLIB = -dynamiclib
DYNLIB_EXT = .dylib
LIB_FLAGS = -install_name @rpath/$(EXE_NAME)$(DYNLIB_EXT)
endif

INC +=  $(INC_UUID) $(INC_SLUG) $(INC_BOOST) $(INC_GSL) $(INC_GRACKLE) $(INC_LAPACKE)
LIB +=  $(LIB_UUID) $(LIB_SLUG) $(LIB_BOOST) $(LIB_GSL) $(LIB_GRACKLE) $(LIB_LAPACKE)

#### General settings after setting system-dependent ones

DEFINES += -DGRACKLE_DIR=\"$(GRACKLE_DIR)\"
DEFINES += -DSLUG_DIR=\"$(SLUG_DIR)\"
LIB_NAME = lib$(EXE_NAME)$(DYNLIB_EXT)

##################################################### 

SRC = $(wildcard *.cxx Initial/*.cxx Hydro/*.cxx Output/*.cxx)
OBJ = $(SRC:.cxx=.o)
DEP = $(SRC:.cxx=.d)
HEADERS = $(wildcard *.H Initial/*.H Hydro/*.H Output/*.H)

############# RULES ########################################
.PHONY: all all_debug exe exe_debug lib lib_debug install copy_headers copy_exe copy_lib clean distclean

all: exe
all: lib

all_debug: exe_debug
all_debug: lib_debug

exe: CXX_FLAGS += $(CXX_FLAGS_OPT)
exe: $(EXE_NAME)

exe_debug: CXX_FLAGS += $(CXX_FLAGS_DEBUG)
exe_debug: $(EXE_NAME)

lib: CXX_FLAGS += $(CXX_FLAGS_OPT)
lib: $(LIB_NAME)

lib_debug: CXX_FLAGS += $(CXX_FLAGS_DEBUG)
lib_debug: $(LIB_NAME)
 
-include $(DEP)

%.o: %.cxx
	$(CXX) $(DEFINES) $(CXX_FLAGS) $(INC) -MMD -MP -c $< -o $@

$(EXE_NAME): $(OBJ)
	$(CXX) $(CXX_FLAGS) -o $(EXE_NAME) $(OBJ) $(LIB)

$(LIB_NAME): $(OBJ)
	$(CXX) $(CXX_FLAGS) -o $(LIB_NAME) $(OBJ) $(LIB) $(DYNLIB) $(LIB_FLAGS)

install: copy_headers copy_exe copy_lib 

copy_headers: 
# structure needs to be preserved manually,
# since BSD and GNU version of `cp` don't have 
# a consistent way of doing it automatically
	@ mkdir -p $(INC_DIR)
	@ mkdir -p $(INC_DIR)/Initial
	@ mkdir -p $(INC_DIR)/Hydro
	@ mkdir -p $(INC_DIR)/Output
	@ mkdir -p $(INC_DIR)/Initial
	@ echo "Copying headers into " $(INC_DIR)
	@ for HEADER in $(HEADERS); \
	do \
		cp -f $$HEADER $(INC_DIR)/$$HEADER ; \
	done

copy_exe: $(EXE_NAME)
	@ mkdir -p $(INC_DIR)
	cp $(EXE_NAME) $(BIN_DIR)

copy_lib: $(LIB_NAME)
	@ mkdir -p $(LIB_DIR)
	cp $(LIB_NAME) $(LIB_DIR)

clean:
	@ rm -f $(OBJ) $(DEP)
	@ rm -f $(EXE_NAME) $(LIB_NAME)
	@ rm -rf *.dSYM

distclean: clean
	@ rm -rf $(INC_DIR)/* 
	@ rm -f  $(BIN_DIR)/*
	@ rm -f  $(LIB_DIR)/*
