# INTEL compiler.
CC = icc

# INTEL compiler flags for OpenMP and optimization.
CFLAGS = -qopenmp -O3 -Wall

# Libraries.
LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Source and binaries directory.
SRC_DIR := ./src
OBJ_DIR := ./bin

# Create build directory if necessary.
$(shell mkdir -p $(OBJ_DIR))

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
	rm -rf $(MAIN) bin
