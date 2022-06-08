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

##
## % 
## GDRUN_TARGETS			: Targets to be run for the grid deformation
GDRUN_TARGETS = Modules Auxiliary IO_output Drivers

##
## % Compiler
## %=========
## FC			: Compiler to be used
FC = gfortran
## CFLAGS			: Compiler flags
CFLAGS = -c -g -Og -Wall

##
## % Files
## %======
## MAIN_RUNFILE			: Main runfile (.F90) - single file
MAIN_RUNFILE = MainRunFileGridDeformation.F90

## DRIVER_FILES			: Driver filenames (.F90) - unsequenced
DRIVER_FILES = $(wildcard src/Drivers/*.F90)

## MODULE_FILES			: Module filenames (.F90) - sequence matters
MODULE_FILES = src/Modules/gdmod_types.F90\
 src/Modules/gdmod_userinput.F90

## AUXILIARY_FILES			: Auxiliary filenames (.F90) - unsequenced
AUXILIARY_FILES = $(wildcard src/Auxiliary/*.F90)

## OUTPUT_FILES			: output generation filenames (.F90) - unsequenced
OUTPUT_FILES = $(wildcard src/IO/Output/*.F90)

## 
