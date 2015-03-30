CC=gcc

FLAGS 		= -O2 -lm 
FLAGS_DEBUG = -ggdb -Wall -lm
DYLD_LIBRARY_PATH = /Users/egentry/bin/grackle/lib  # ADD LINUX FRIENDLY OPTION!

EXE_FILE 	= fluid.adiabatic

OBJS 		= fluid.adiabatic.o

############# SUFFIX RULES ########################################

$(EXE_FILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(FLAGS) $(OBJS) -o $(EXE_FILE)

.SUFFIXES :             # clear all defaults
.SUFFIXES : .c .o       # and replace

.c.o:
	$(CC) $(FLAGS) -c $<

############# MAKE RULES ########################################

fluid.adiabatic: fluid.adiabatic.c
	$(CC) -O2 -c fluid.adiabatic.c -lm
	$(CC) -O2 -o fluid.adiabatic fluid.adiabatic.o -lm

clean:
	@rm -f *.o
	@rm -f $(EXE_FILE)


debug: FLAGS=$(FLAGS_DEBUG)
debug: clean     				# is there a better way to force .o files to be gdb friendly?
debug: $(EXE_FILE)


############# DEPENDENCIES ########################################
