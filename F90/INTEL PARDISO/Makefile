# INTEL compiler.
CC = icc
F90 = ifort

# INTEL compiler flags for OpenMP and optimization.
CFLAGS = -qopenmp -O3 -Wall
F90FLAGS = -qopenmp -O3

# Libraries.
LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lstdc++

# Source and binaries directory.
SRC_DIR := ./src
OBJ_DIR := ./bin

# Create build directory if necessary.
$(shell mkdir -p $(OBJ_DIR))

# Program sources and binaries.
C_SRCS := $(wildcard $(SRC_DIR)/*.cpp)
F_SRCS := $(wildcard $(SRC_DIR)/*.f90)
SRCS := $(C_SRCS) $(F_SRCS)
C_OBJS := $(subst src/,bin/,$(subst .cpp,.o,$(C_SRCS))) 
F_OBJS := $(subst src/,bin/,$(subst .f90,.o,$(F_SRCS))) 
OBJS := $(C_OBJS) $(F_OBJS)

# Main execuatble object dependencies.
MAIN = ELLSOLVERF

all: $(MAIN)
	@echo Executable has been compiled succesfully.

bin/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
bin/%.o: src/%.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(MAIN): $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) -o $(MAIN) $(LIBS)
clean:
	rm -rf $(MAIN) bin
