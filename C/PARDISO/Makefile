# INTEL compiler.
CC = icc

# INTEL compiler flags for OpenMP and optimization.
CFLAGS = -qopenmp -O3 -Wall

# Libraries.
LIBS = -lpardiso500-INTEL1301-X86-64 -lblas -llapack -lmkl_intel_lp64 -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -liomp5 -lgfortran -fopenmp -lpthread -lm -lstdc++ -lifcore

# Source and binaries directory.
SRC_DIR := ./src
OBJ_DIR := ./bin

# Program sources and binaries.
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(subst src/,bin/,$(subst .cpp,.o,$(SRCS))) 

# Main execuatble object dependencies.
MAIN = ELLSOLVERC

all: $(MAIN)
	@echo Executable has been compiled succesfully.

bin/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(MAIN) $(LIBS)
clean:
	rm -f $(MAIN) bin/*.o
