# Configuration file for the grid deformation makefile. Variables etc 
# are listed below. When introducing new variables or other 
# functionality, please document using the '##' prefix. THis should
# be compatible with the 'help' target in the makefile. 

##
## %===================================================================%
## %                                                                   %
## %                             VARIABLES                             %
## %                                                                   %
## %===================================================================%
##
## Variables are defined in config.mk

## % Library paths
## %==============
## LAPACKPATH 			: LAPACK library path (user defined)
LAPACKPATH = -lopenblas

## BLASPATH 			: BLAS library path (user defined)
BLASPATH = -lopenblas

## UMFPACKPATH 			: UMFPACK library path (user defined)
UMFPACKPATH = -lumfpack

## % Include paths
## %==============
## SUITESPARSEPATH      : SuiteSparse header file path
SUITESPARSEPATH = /usr/include/suitesparse

##
## % Targets
## %========
## GDRUN_TARGETS			: Targets to be run for the grid deformation
GDRUN_TARGETS = Clayer ClayerF Constants Auxiliary General Numerics Optimization Modules IO_b25  \
    IO_carre IO_output IO_input  Setup  Drivers 

## GOAT_TARGETS             : Targets to be run for the full goat
GOAT_TARGETS = $(GDRUN_TARGETS) 

## TEST_TARGETS             : Targets to be run for goat tests
TEST_TARGETS = $(GOAT_TARGETS) 

## CTEST_TARGETS            : Targets to be run to test C layer
CTEST_TARGETS = Clayer

##
## % Compiler
## %=========
## FC			: Compiler to be used for fortran
FC = gfortran
## CFLAGS			: Compiler flags for standard compilation (may be overridden)
CFLAGS_DEF = -c -g -Wall -O0 -Wno-unused-dummy-argument -Wno-maybe-uninitialized
## CFLAGS_OMP	: compiler flags for OpenMP 
CFLAGS_OMP = -c -Wall -fopenmp
## CFLAGS_DEBUG		: compiler flags for debugging
CFLAGS_DEBUG = -c -g -Wall -pg  -O0 
## CFLAGS_OMP_DEBUG	: compiler flags for OpenMP and debugging
CFLAGS_OMP_DEBUG = -c -g -Wall -pg  -O0 -fopenmp

## CC           : Compiler to be used for C
CC = gcc 

CCFLAGS_DEF = -c -g -Wall -O0

## % Linker
## %=======
## LFLAGS			: linking flags to be used (apart from libraries)
LFLAGS_DEF =  

## LFLAGS_DEBUG 	: linking flags for debugging
LFLAGS_DEBUG = -pg -g

## LFLAGS			: linking flags to be used for openMP
LFLAGS_OMP = -pg -g -fopenmp

## LFLAGS_DEBUG 	: linking flags for debugging
LFLAGS_OMP_DEBUG = -pg -g -fopenmp

# Set flags
#==========
# Set the CFLAGS
CFLAGS = $(CFLAGS_DEF)

# Set the linking flags
LFLAGS = $(LFLAGS_DEF)

# Set CFLAGS for C compiler
CCFLAGS = $(CCFLAGS_DEF) 

##
## % Files
## %======
## MAIN_RUNFILE			: Main runfile (.F90) - single file
MAIN_RUNFILE = MainRunFileGridDeformation.F90

## GENERAL_FILES				: All general files (e.g. precision definition, ... )
GENERAL_FILES = src/General/mod_sparseinterface.F90 src/General/mod_readwrite.F90 $(wildcard src/General/*.F90)
    

## DRIVER_FILES			: Driver filenames (.F90) - unsequenced
DRIVER_FILES = $(wildcard src/Drivers/*.F90)

## MODULE_FILES			: Module filenames (.F90, .F) - sequence matters
MODULE_FILES = $(wildcard src/Modules/Goat/*.F90)\
    src/Modules/GD/gdmod_types.F90 src/Modules/GD/gdmod_userinput.F90 src/Modules/GD/gdmod_plots.F90 src/Modules/GD/gdmod_designvariables.F90 \
    src/Modules/GD/gdmod_utility_optimization.F90 src/Modules/GD/gdmod_constraints.F90\
    $(wildcard src/Modules/GD/*.F90) \
    $(wildcard src/Modules/*.F90) $(wildcard src/Modules/*.F)


## AUXILIARY_FILES			: Auxiliary filenames (.F90) - unsequenced
AUXILIARY_FILES = src/Auxiliary/mod_plotter.F90 \
    src/Auxiliary/Construct2DStructuredGrid.F90 \
    $(wildcard src/Auxiliary/*.F90) \
    src/Auxiliary/Interpolation/Interpolant2D_auxiliaries.F90 \
    src/Auxiliary/Interpolation/Interpolant2D.F90 \
    src/Auxiliary/Interpolation/BicubicSplineInterpolant.F90 \
    src/Auxiliary/Interpolation/StructuredInterpolant2D.F90 \
    $(wildcard src/Auxiliary/Interpolation/*.F90) 

## B25_FILES			: b25 generation filenames (.F90, .F) - unsequenced
B25_FILES = $(wildcard src/IO/B25/*.F90) $(wildcard src/IO/B25/*.F)

## CARRE_FILES			: carre generation filenames (.F90, .F) - unsequenced
CARRE_FILES = $(wildcard src/IO/CARRE/*.F90) $(wildcard src/IO/CARRE/*.F)

## INPUT_FILES			: input generation filenames (.F90) - unsequenced
INPUT_FILES = $(wildcard src/IO/Input/*.F90)

## OUTPUT_FILES			: output generation filenames (.F90) - unsequenced
OUTPUT_FILES = $(wildcard src/IO/Output/*.F90)

## SETUP_FILES			: setup file generation names (.F90) - unsequenced
SETUP_FILES = $(wildcard src/Setup/*.F90)

## NUMERICS_FILES 		: numeric file generation names (.F90) - unsequenced
NUMERICS_FILES = src/Numerics/PolygonLevelsetFunction2D.F90 $(wildcard src/Numerics/*.F90)

## OPTIMIZATION 		: optimization file generation names (.F90) - unsequenced
OPTIMIZATION_FILES = src/Optimization/optmod_designvariables.F90 src/Optimization/optmod_costfunction.F90 \
    src/Optimization/optmod_constraints.F90 src/Optimization/optmod_monitor.F90 \
    src/Optimization/optmod_numerics.F90 $(wildcard src/Optimization/*.F90)

## CONSTANTS            : constants such as precision and special characters (.F90) - unsequenced
CONSTANTS_FILES = src/Constants/mod_precision.F90 $(wildcard src/Constants/*.F90)

## Clayer               : c files for interfacing with other c code
CLAYER_FILES    = $(wildcard src/Clayer/*.c)

## ClayerF              : fortran files for interfacing with other c code
CLAYERF_FILES    =  src/Clayer/CSparseF.F90 src/Clayer/Clayer.F90
