# Include solver sources.


# Help or wrong option print-out.
help:
	@echo $$'$(WRONG_OPTION)'
	@echo "Elliptic Solver compile help."
	@echo ""
	@echo "Usage: make Target [Options...]"
	@echo ""
	@echo "   Target:"
	@echo "      C		- Solver main program and subroutine calls are C based."
	@echo "      FORTRAN	- Solver main program and subroutine calls are FORTRAN based."
	@echo "      help	- Print this help."
	@echo ""
	@echo "   Options:"
	@echo "      compiler={intel|gnu}"
	@echo "          Specifies whether to use Intel compilers or GNU's alternatives."
	@echo "          Default: intel."
	@echo "      pardiso={mkl|vanilla}"
	@echo "          Specifies whether to use MKL's version of PARDISO or the official one."
	@echo "          Default: mkl."
#	@echo "      MKLROOT=<MKL_directory>"
#	@echo "          Specifies the location of MKL libraries used."
#	@echo "          Default: The Intel MKL installation directory."
	@echo ""
	@echo "Usage examples:"
	@echo ""
	@echo "   make C compiler=gnu pardiso=vanilla"
	@echo "      Create C based executable using gcc and PARDISO's official version."
	@echo ""
	@echo "   make FORTRAN compiler=intel pardiso=mkl"
	@echo "      Create FORTRAN based executable using ifort, icc and MKL's PARDISO."
	@echo "" 
# No MKL installation.
#no_mkl:
#	@echo ""
#	@echo ""
#	@echo "*** CRITICAL ERROR: There does not seem to be an MKL installation."
#	@echo ""
#	@echo ""
#	@echo "Please install MKL and configure it properly with the scripts provided in"
#	@echo "the installlation. A succesfull installation will have MKLROOT in your"
#	@echo "enviroment variables."
#	@echo ""

# Defaults
compiler=intel
interface=lp64
threading=omp
pardiso=mkl

# Check for compiler option.
ifneq ($(compiler),intel)
ifneq ($(compiler),gnu)
MSG2+=compiler=$(compiler)
endif
endif

# Check for PARDISO version option.
ifneq ($(pardiso),mkl)
ifneq ($(pardiso),vanilla)
MSG2+=pardiso=$(pardiso)
endif
endif

# Check for errors in command line options.
ifneq ("$(MSG2)","")
WRONG_OPTION=\n\n*** COMMAND LINE ERROR: Wrong value of option(s): $(MSG2)\n\n
TARGET=help
endif

# Actual compiling.
ifdef _IA

# Static vs. Dynamic linking.
ifdef ($(SD),static)
 EXT=a
 RES_EXT=lib
 S_GRP=-Wl,--start-group
 E_GRP=-WL,--end-group
else
 EXT=so
 RES_EXT=so
 S_GRP=
 E_GRP=
endif

# Check for MKLROOT
ifndef MKLROOT
 $(warning *** MKLROOT is not setup. Perhaps user has not installed MKL or set it up properly.)
 $(error Try >make help)
endif
MKL_PATH = $(MKLROOT)/lib/$(_IA)
CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(_IA)
TBB_PATH=$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.4

# Linker options.
LOPTS = 

# Compiler options.
ifeq ($(compiler),gnu)
 override COMPILER=gcc
 IFACE_COMP_PART=intel
 IFACE_THREADING_PART=gnu
 ifeq ($(RES_EXT),so)
  LOPTS = -Wl,--no-as-needed
 endif
else
 override COMPILER=icc
 IFACE_COMP_PART=intel
 IFACE_THREADING_PART=intel
endif
OPTIONS = -w

# PARDISO options.
ifeq ($(pardiso),vanilla)
 PARDISO_OPTS = 

# Architecture options.
ifeq ($(_IA),ia32)
 ifeq ($(interface),ilp64)
  $(warning *** ILP64 interface is not available for MKL-IA32)
  $(error Try >make help)
 endif
 if ($(compiler),intel)
  # This option tells the compiler to generate optimized code for Pentium or later processor.
  OPTIONS += -mia32
 endif
 IFACE_SUFF=
 M32_64 = -m64 # This option tells the compiler to generate code for IA-32 architecture.
else
 IFACE_SUFF=_$(threading)
 M32_64=-m64 # This option tells the compiler to generate code for Intel64 architecture.
endif

# MKL library paths.
ifeq ($(EXT),so)
 MKL_LD_PATH=-L"$(MKL_PATH)"
 MKL_PREFIX=-l
 MKL_SUFFIX=
else
 MKL_LD_PATH=
 MKL_PREFIX="$(MKL_PATH)/lib
 MKL_SUFFIX=.a"
endif
IFACE_LIB=$(MKL_PREFIX)mkl_$(IFACE_COMP_PART)$(IFACE_SUFF)$(MKL_SUFFIX)

# Interface options.
ifeq ($(interface),ilp64)
 OPTIONS += -DMKL_ILP64
endif

# Threading options.
ifeq ($(threading),sequential)
 THREADING_LIB=$(MKL_PREFIX)mkl_sequential$(MKL_SUFFIX)
 OMP_LIB = 
else
ifeq ($(threading),omp)
 THREADING_LIB=$(MKL_PREFIX)mkl_$(IFACE_THREADING_PART)_thread$(MKL_SUFFIX)
 OMP_LIB = -L"$(CMPLR_PATH)" -liomp5
else
 THREADING_LIB=$(MKL_PREFIX)mkl_tbb_thread$(MKL_SUFFIX)
 OMP_LIB = -L"$(TBB_PATH)" -ltbb -lstdc++
endif
endif

CORE_LIB=$(MKL_PREFIX)mkl_core$(MKL_SUFFIX)

endif # ifdef _IA

#------------------------------------------------------------------------------
# C Target compilation.
#------------------------------------------------------------------------------
C: 
	@$(MAKE) $(TARGET) --no-print-directory SD=dynamic _IA=intel64
