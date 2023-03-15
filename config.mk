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
LAPACKPATH = /data/leuven/320/vsc32061/lapack/lapack/liblapack.a
## BLASPATH 			: BLAS library path (user defined)
BLASPATH = /data/leuven/320/vsc32061/lapack/lapack/librefblas.a
## UMFPACKPATH 			: UMFPACK library path (user defined)
UMFPACKPATH = -lumfpack

##
## % Targets
## %========
## GDRUN_TARGETS			: Targets to be run for the grid deformation
GDRUN_TARGETS = General Optimization Modules Auxiliary IO_b25 IO_carre \
    IO_output IO_input  Setup  Drivers Numerics

## GOAT_TARGETS             : Targets to be run for the full goat
GOAT_TARGETS = $(GDRUN_TARGETS) 

##
## % Compiler
## %=========
## FC			: Compiler to be used
FC = gfortran
## CFLAGS			: Compiler flags for standard compilation (may be overridden)
CFLAGS = -c -Wall -O3
## CFLAGS_OMP	: compiler flags for OpenMP 
CFLAGS_OMP = -c -Wall -fopenmp
## CFLAGS_DEBUG		: compiler flags for debugging
CFLAGS_DEBUG = -c -g -Wall -pg  -O0 
## CFLAGS_OMP_DEBUG	: compiler flags for OpenMP and debugging
CFLAGS_OMP_DEBUG = -c -g -Wall -pg  -O0 -fopenmp

## % Linker
## %=======
## LFLAGS			: linking flags to be used (apart from libraries)
LFLAGS =  

## LFLAGS_DEBUG 	: linking flags for debugging
LFLAGS_DEBUG = -pg -g

## LFLAGS			: linking flags to be used for openMP
LFLAGS_OMP = -pg -g -fopenmp

## LFLAGS_DEBUG 	: linking flags for debugging
LFLAGS_OMP_DEBUG = -pg -g -fopenmp

# Set flags
#==========
# Set the CFLAGS
CFLAGS = $(CFLAGS_DEBUG)

# Set the linking flags
LFLAGS = $(LFLAGS_DEBUG)

##
## % Files
## %======
## MAIN_RUNFILE			: Main runfile (.F90) - single file
MAIN_RUNFILE = MainRunFileGridDeformation.F90

## GENERAL_FILES				: All general files (e.g. precision definition, ... )
GENERAL_FILES = $(wildcard src/General/*.F90)

## DRIVER_FILES			: Driver filenames (.F90) - unsequenced
DRIVER_FILES = $(wildcard src/Drivers/*.F90)

## MODULE_FILES			: Module filenames (.F90, .F) - sequence matters
MODULE_FILES =  $(wildcard src/Modules/*.F90) $(wildcard src/Modules/*.F)

## AUXILIARY_FILES			: Auxiliary filenames (.F90) - unsequenced
AUXILIARY_FILES = $(wildcard src/Auxiliary/*.F90)

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
NUMERICS_FILES = $(wildcard src/Numerics/*.F90)

## OPTIMIZATION 		: optimization file generation names (.F90) - unsequenced
OPTIMIZATION_FILES = $(wildcard src/Optimization/*.F90)
