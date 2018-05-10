# -----------------------------------------------------------------------------
# DEFAULTS AND INITIALIZATION.
# -----------------------------------------------------------------------------
# C compiler.
CC = icc
# FORTRAN compiler.
F90 = ifort

# C flags for 64 bit architecture and optimization.
CFLAGS = -m64 -O2
# FORTRAN flags: same.
F90FLAGS = -m64 -O2

# Linker flags.
LDFLAGS = 

# Preprocessor.
FORTRAN_PP = 

# MKL: Default are INTEL setups.
# MKL libraries path.
MKL_PATH=$(MKLROOT)/lib/intel64
# Linker instruction path.
MKL_LD_PATH=-L"$(MKL_PATH)"
# MKL libraries.
MKL_MAIN_LIBS = -lmkl_core -lmkl_intel_thread
MKL_INTEL_LIB = -lmkl_intel_lp64
MKL_FORTRAN_LIB  = 

# OpenMP libraries.
OMP_LIBS = -liomp5

# Other libraries.
OTHER_LIBS = -lpthread -lm -ldl
FORTRAN_LIBS = -lstdc++

# -----------------------------------------------------------------------------
# HELP AND SANITY CHECKS.
# -----------------------------------------------------------------------------
help:
	@echo "Elliptic Solver help."
	@echo ""
	@echo "Usage make Target [Options...]"
	@echo ""
	@echo "   Target:"
	@echo "      C          - Compile solver using C main program."
	@echo "      FORTRAN    - Compile solver using FORTRAN main program."
	@echo "      clean      - Remove binaries and executable."
	@echo "      help       - Print this help."
	@echo ""
	@echo "   Options:"
	@echo "      compiler={intel|gnu}"
	@echo "         Specifies whether to use Intel's icc or GNU's gcc."
	@echo "         Default: intel."
	@echo ""


# Check compiler options.
ifneq ($(compiler),intel)
  ifneq ($(compiler),gnu)
    MSG += compiler = $(compiler)
  endif
endif

# Check for errors in command line options.
ifneq ("$(MSG)","")
  WRONG_OPTION = \n\n*** COMMAND LINE ERROR: Wrong value of option(s): $(MSG)\n\n
  TARGET = help
endif

# -----------------------------------------------------------------------------
# SETUP VARIABLES.
# -----------------------------------------------------------------------------
# Set compiler first.
ifeq ($(compiler),gnu)
  override CC = gcc
  override F90 = gfortran
else
  ifeq ($(compiler),intel)
    override CC = icc
    override F90 = ifort
  endif
endif

# Setup options.
ifeq ($(compiler),gnu)
  # Modify flags for OpenMP.
  CFLAGS += -fopenmp
  F90FLAGS += -fopenmp
  # Modify linker flags.
  LDFLAGS += -Wl,--no-as-needed
else
  # Intel OpenMP.
  CFLAGS += -qopenmp
  F90FLAGS += -qopenmp
endif

# Check for gfortran compiler.
ifeq ($(compiler),gnu)
  MKL_FORTRAN_LIB = -lmkl_gf_lp64
else
  MKL_FORTRAN_LIB = -lmkl_intel_lp64
endif

# -----------------------------------------------------------------------------
# SOURCES AND BINARIES.
# -----------------------------------------------------------------------------
# Source and binaries directory.
SRC_DIR := ./src
OBJ_DIR := ./bin

# Create build directory if necessary.
$(shell mkdir -p $(OBJ_DIR))

# Program sources and binaries.
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(subst src/,bin/,$(subst .cpp,.o,$(SRCS))) 
C_MAIN_SRC := src/main.cpp
F_MAIN_SRC := src/main.f90
C_SRCS := src/elliptic_tools.cpp src/flat_laplacian.cpp src/flat_laplacian_csr_gen.cpp src/general_elliptic.cpp src/general_elliptic_csr_gen.cpp src/low_rank.cpp src/pardiso_start.cpp src/pardiso_stop.cpp src/pardiso_wrapper.cpp src/tools.cpp

C_MAIN_OBJ := bin/main_c.o
F_MAIN_OBJ := bin/main_f.o
C_OBJS := bin/elliptic_tools.o bin/flat_laplacian.o bin/flat_laplacian_csr_gen.o bin/general_elliptic.o bin/general_elliptic_csr_gen.o bin/low_rank.o bin/pardiso_start.o bin/pardiso_stop.o bin/pardiso_wrapper.o bin/tools.o

# -----------------------------------------------------------------------------
# MAIN COMPILATION AND LINKING.
# -----------------------------------------------------------------------------
#  Executable name.
C_EXE = ELLSOLVEC
F_EXE = ELLSOLVEF

# C-based executable.
C: $(C_EXE)

# FORTRAN-based executable.
FORTRAN: FORTRAN_PP = -D FORTRAN
FORTRAN: $(F_EXE)

# C main file.
$(C_MAIN_OBJ): $(C_MAIN_SRC)
	@echo ""
	@echo "Compiling C main program..."
	$(CC) $(CFLAGS) -c $< -o $@

# FORTRAN main file.
$(F_MAIN_OBJ): $(F_MAIN_SRC)
	@echo ""
	@echo "Compiling FORTRAN main program..."
	$(F90) $(F90FLAGS) -c $< -o $@

# Program C binaries.
bin/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(FORTRAN_PP) -c $< -o $@

# Link C executable.
$(C_EXE): $(C_MAIN_OBJ) $(C_OBJS)	
	@echo ""
	@echo "Linking with C compiler..."
	$(CC) $(CFLAGS) $(C_OBJS) $(C_MAIN_OBJ) -o $(C_EXE) $(LDFLAGS) $(MKL_LD_PATH) $(MKL_MAIN_LIBS) $(MKL_INTEL_LIB) $(OMP_LIBS) $(OTHER_LIBS)

# Link FORTRAN executable.
$(F_EXE): $(F_MAIN_OBJ) $(C_OBJS)
	@echo ""
	@echo "Linking with FORTRAN compiler..."
	$(F90) $(F90FLAGS) $(C_OBJS) $(F_MAIN_OBJ) -o $(F_EXE) $(LDFLAGS) $(MKL_LD_PATH) $(MKL_MAIN_LIBS) $(MKL_FORTRAN_LIB) $(OMP_LIBS) $(OTHER_LIBS) $(FORTRAN_LIBS)

# Clean up binaries and executable.
clean:
	@echo "Cleaning up executables and binaries..."
	rm -rf $(C_EXE) $(F_EXE) bin
