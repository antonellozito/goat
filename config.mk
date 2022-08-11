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
GDRUN_TARGETS = Modules Auxiliary IO_b25 IO_output IO_input  Setup  Drivers Numerics

##
## % Compiler
## %=========
## FC			: Compiler to be used
FC = gfortran
## CFLAGS			: Compiler flags
CFLAGS = -c -g -Og -Wall -pg 

##
## % Files
## %======
## MAIN_RUNFILE			: Main runfile (.F90) - single file
MAIN_RUNFILE = MainRunFileGridDeformation.F90

## DRIVER_FILES			: Driver filenames (.F90) - unsequenced
DRIVER_FILES = $(wildcard src/Drivers/*.F90)

## MODULE_FILES			: Module filenames (.F90, .F) - sequence matters
MODULE_FILES = $(wildcard src/Modules/*.F90) $(wildcard src/Modules/*.F)

## AUXILIARY_FILES			: Auxiliary filenames (.F90) - unsequenced
AUXILIARY_FILES = $(wildcard src/Auxiliary/*.F90)

## B25_FILES			: b25 generation filenames (.F90, .F) - unsequenced
B25_FILES = $(wildcard src/IO/B25/*.F90) $(wildcard src/IO/B25/*.F)

## INPUT_FILES			: input generation filenames (.F90) - unsequenced
INPUT_FILES = $(wildcard src/IO/Input/*.F90)

## OUTPUT_FILES			: output generation filenames (.F90) - unsequenced
OUTPUT_FILES = $(wildcard src/IO/Output/*.F90)

## SETUP_FILES			: setup file generation names (.F90) - unsequenced
SETUP_FILES = $(wildcard src/Setup/*.F90)

## NUMERICS_FILES 		: numeric file generation names (.F90) - unsequenced
NUMERICS_FILES = $(wildcard src/Numerics/*.F90)
