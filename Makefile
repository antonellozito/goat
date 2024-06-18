# A very simple and straightforward makefile without too many options 
# yet. To be improved in the future. 

# Include the config file
include config.mk

##
## %===================================================================%
## %                                                                   %
## %                             TARGETS                               %
## %                                                                   %
## %===================================================================%
##
## % Main targets
## %=============
# Linking
## goat			: Create executable for goat
## tests		: Create executable for tests
## testc 		: Create executable for tests of C layer
## goattranslator: Create executable for GOAToptions file translator
## shapeopt 	: Create executable for shape optimization with goat

goat: $(GOAT_TARGETS) goat.o
	$(FC) -o goat *.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include

goattranslator: $(GOATTRANSLATOR_TARGETS) goattranslator.o 
	$(FC) -o goattranslator *.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include

tests: $(TEST_TARGETS) tests.o 
	$(FC) -o tests *.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include

gdrun: $(GDRUN_TARGETS) MainRunFileGridDeformation.o
	$(FC) -o gdrun *.o $(LFLAGS) -l $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH)

testc: $(CTEST_TARGETS) testc.o 
	$(CC) -o testc *.o -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include

shapeopt: $(SHAPEOPT_TARGETS) shapeopt.o 
	$(FC) -o shapeopt *.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include

shapeopt_solps: $(SHAPEOPTSOLPS_TARGETS) shapeopt.o 
	$(FC) -o shapeopt *.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include -I $(B25LIBPATH) -I $(B25ADJLIBPATH) $(DEFINEFLAGS)


## % Runfiles
## %=========
## MainRunFileGridDeformation.o			: main runfile grid deformation
# Compiling
MainRunFileGridDeformation.o: Runfiles/MainRunFileGridDeformation.F90
	$(FC) $(CFLAGS) Runfiles/MainRunFileGridDeformation.F90 

## Goat.o 			: main run file goat
goat.o: Runfiles/Goat.F90
	$(FC) $(CFLAGS) Runfiles/Goat.F90 

## Tests.o 			: all tests
tests.o: Runfiles/Tests.F90 
	$(FC) $(CFLAGS) Runfiles/Tests.F90 

## Testc.o 			: C layer tests
testc.o: Runfiles/Testc.c 
	$(CC) $(CCFLAGS) Runfiles/Testc.c -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include

## goattranslator.o	: translator
goattranslator.o: Runfiles/TranslateGOAToptionsFile.F90
	$(FC) $(CFLAGS) Runfiles/TranslateGOAToptionsFile.F90 

## shapeopt.o		: shape optimization 
shapeopt.o : Runfiles/ShapeOptimization.F90 
	$(FC) $(CFLAGS) Runfiles/ShapeOptimization.F90 $(DEFINEFLAGS)

##
## % Folder compilation targets
## %===========================
## Constants 		: compile files containing constants
Constants: $(CONSTANTS_FILES) 
	$(FC) $(CFLAGS) $^
	touch Constants 

## General 			: compile general files and modules
General: $(GENERAL_FILES)
	$(FC) $(CFLAGS) $^
	touch General
	
## Modules			: compile modules 
Modules: $(MODULE_FILES)
	$(FC) $(CFLAGS) $^
	touch Modules

## Auxiliary			: compile auxiliary routines
Auxiliary: $(AUXILIARY_FILES)
	$(FC) $(CFLAGS) $^
	touch Auxiliary

## Drivers			: compile all driver routines
Drivers: $(DRIVER_FILES)
	$(FC) $(CFLAGS) $^
	touch Drivers

## Setup			: compile all setup routines
Setup: $(SETUP_FILES)
	$(FC) $(CFLAGS) $^
	touch Setup

## IO_output 			: compile output routines
IO_output: $(OUTPUT_FILES)
	$(FC) $(CFLAGS) $^
	touch IO_output

## IO_input			: compile input routines
IO_input: $(INPUT_FILES)
	$(FC) $(CFLAGS) $^
	touch IO_input

## IO_b25			: compile b25 routines
IO_b25: $(B25_FILES)
	$(FC) $(CFLAGS) $^
	touch IO_b25

## IO_carre			: compile carre routines
IO_carre: $(CARRE_FILES)
	$(FC) $(CFLAGS) $^
	touch IO_carre

## Numerics		 	: compile numerics routines
Numerics: $(NUMERICS_FILES)
	$(FC) $(CFLAGS) $^
	touch Numerics

## Optimization		 	: compile optimization routines
Optimization: $(OPTIMIZATION_FILES)
	$(FC) $(CFLAGS) $^
	touch Optimization

## Clayer 				: compile C interlayer routines
Clayer: $(CLAYER_FILES)
	$(CC) $(CCFLAGS) $^ -I $(SUITESPARSEPATH) -I src/Clayer/Include
	touch Clayer

## ClayerF 				: compile fortran modules in Clayer
ClayerF: $(CLAYERF_FILES)
	$(FC) $(CFLAGS) $^
	touch ClayerF

## ShapeOptimization 			: compile shape optimization modules
ShapeOptimization: $(SHAPEOPT_FILES)
	$(FC) $(CFLAGS) $(DEFINEFLAGS) $^
	touch ShapeOptimization

## ShapeOptimizationSolps 			: compile shape optimization modules for SOLPS
ShapeOptimizationSolps: $(SHAPEOPTSOLPS_FILES)
	$(FC) $(CFLAGS) -I$(B25LIBPATH) -I$(B25ADJLIBPATH) $(DEFINEFLAGS) $^ 
	touch ShapeOptimizationSolps


##
## % Run commands
## %=============
## gd			: Run grid deformation with gdrun
# Run
.PHONY: gd
gd: 
	gdrun
	make clean

##
## % Auxiliary targets
## %==================
## clean			: clean by removing *.o 
# Cleanup
.PHONY: clean
clean: 
	rm *.o $(wildcard gdrun*); rm $(GDRUN_TARGETS); rm $(GOATTRANSLATOR_TARGETS); rm $(SHAPEOPT_TARGETS); rm $(SHAPEOPTSOLPS_TARGETS);

## deepclean			: clean by removing *.o, *.mod, and executables
# Cleanup
.PHONY: deepclean
deepclean: 
	rm *.o *.mod $(wildcard gdrun*); \
	find . -name "*.mod" -type f -delete; \
	rm $(GOAT_TARGETS); \
	rm $(GOATTRANSLATOR_TARGETS); \
	rm $(SHAPEOPT_TARGETS); \
	rm $(SHAPEOPTSOLPS_TARGETS); \
	rm goat; \
	rm tests; \
	rm goattranslator; \
	rm testc; \
	rm shapeopt; \

## cleanexeo			: clean o-files of executables
.PHONY: cleanexeo
cleanexeo:	
	rm Goat.o; rm Tests.o; rm Testc.o; rm TranslateGOAToptionsFile.o; rm shapeopt.o
	

## help			: print out documentation
# Help - prints out all the ## statements
.PHONY : help 
help : config.mk Makefile  
	@sed -n 's/^##//p' $^

##
##