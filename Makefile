# A very simple and straightforward makefile without too many options 
# yet. To be improved in the future. An attempt has been made to clean
# up the compilation and linking steps etc by moving all .o and .mod
# files into a build directory that is constructed during compilation

# Include the config file
include config.mk

# Echo the different environment variables 
$(info % ===========================)
$(info % Executing the goat makefile)
$(info % ===========================)
$(info % )
$(info % General environment variables)
$(info % -----------------------------)
$(info % fortran compiler: $(FC))
$(info % C compiler: $(CC))
$(info % Build directory: $(BUILDDIR))
$(info % )
$(info % SOLPS specific environment variables)
$(info % ------------------------------------)
$(info % running at host: $(HOST_NAME))
$(info % solps path: $(SOLPSTOP))
$(info % B2.5 build path: $(BUILDDIR))
$(info % B2.5 libary path: $(B25LIBPATH))

## %===================================================================%
## %                                                                   %
## %                          PRELIMINARIES                            %
## %                                                                   %
## %===================================================================%
## Construction of directories
## %==========================
## CREATE_BUILDDIR: construct build directory based on BUILDDIR environment variable
CREATE_BUILDDIR:= $(shell mkdir ./builds)
CREATE_BUILDDIR:= $(shell mkdir ./builds/$(BUILDDIR))

## CREATE_EXEDIR: construct executable directory 
CREATE_EXEDIR:= $(shell mkdir ./executables)

BUILDDIR :=./builds/$(BUILDDIR)


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

goat: $(addprefix $(BUILDDIR)/, $(GOAT_TARGETS)) $(BUILDDIR)/goat.o
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(FC) $(LFLAGS) -o $(BUILDDIR)/goat.exe $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include; 
	rm $(BUILDDIR)/Goat.o; 
	cp $(BUILDDIR)/goat.exe ./executables/.

goattranslator: $(addprefix $(BUILDDIR)/,$(GOATTRANSLATOR_TARGETS) ) $(BUILDDIR)/goattranslator.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(FC) $(LFLAGS) -o $(BUILDDIR)/goattranslator.exe $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include
	rm $(BUILDDIR)/TranslateGOAToptionsFile.o; 
	cp $(BUILDDIR)/goattranslator.exe ./executables/.

tests: $(addprefix $(BUILDDIR)/,$(TEST_TARGETS) ) $(BUILDDIR)/tests.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(FC) $(LFLAGS) -o $(BUILDDIR)/tests.exe $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include
	rm $(BUILDDIR)/tests.o; 
	cp $(BUILDDIR)/tests.exe ./executables/.

gdrun: $(addprefix $(BUILDDIR)/,$(GDRUN_TARGETS) ) $(BUILDDIR)/MainRunFileGridDeformation.o
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(FC) $(LFLAGS) -o $(BUILDDIR)/gdrun.exe $(BUILDDIR)/*.o $(LFLAGS) -l $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH)
	rm $(BUILDDIR)/gdrun.o; 
	cp $(BUILDDIR)/gdrun.exe ./executables/.

testc: $(addprefix $(BUILDDIR)/,$(CTEST_TARGETS) ) $(BUILDDIR)/testc.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(CC) $(LFLAGS) -o $(BUILDDIR)/testc.exe $(BUILDDIR)/*.o -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include
	rm $(BUILDDIR)/Testc.o; 
	cp $(BUILDDIR)/testc.exe ./executables/.

shapeopt: $(addprefix $(BUILDDIR)/, $(SHAPEOPT_TARGETS) ) $(BUILDDIR)/shapeopt.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
ifdef DOSOLPS
	$(FC) $(LFLAGS) -o $(BUILDDIR)/shapeopt.exe $(BUILDDIR)/*.o $(B25LIBBPATH)/adStack.o \
	 $(B25LIBBPATH)/b2mod_cdf.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) \
	 -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include  -I$(B25LIBPATH) -L$(B25LIBPATH) -l:libb2.a -L$(B25LIBPATH) -l:libb2.a -lnetcdf $(LD_NETCDF)
else
	$(FC) $(LFLAGS) -o $(BUILDDIR)/shapeopt.exe $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) -lcxsparse \
	-I $(SUITESPARSEPATH) -I src/Clayer/Include
endif
	rm $(BUILDDIR)/ShapeOptimization.o; 
	cp $(BUILDDIR)/shapeopt.exe ./executables/.

#shapeopt_solps: $(addprefix $(BUILDDIR)/,$(SHAPEOPTSOLPS_TARGETS) ) $(BUILDDIR)/shapeopt_solps.o 
#	-mv -f *.o $(BUILDDIR);  
#	-mv -f *.mod $(BUILDDIR); 
#	$(FC) $(LFLAGS) -o $(BUILDDIR)/shapeopt_solps $(BUILDDIR)/*.o $(B25LIBBPATH)/adStack.o \
	 $(B25LIBBPATH)/b2mod_cdf.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) \
	 -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include  -I$(B25LIBPATH) -L$(B25LIBPATH) -l:libb2.a -L$(B25LIBPATH) -l:libb2.a -lnetcdf $(LD_NETCDF)
#	rm $(BUILDDIR)/ShapeOptimization.o; 
#	cp $(BUILDDIR)/shapeopt_solps.exe ./executables/.


## % Runfiles
## %=========
## MainRunFileGridDeformation.o			: main runfile grid deformation
# Compiling
$(BUILDDIR)/MainRunFileGridDeformation.o: Runfiles/MainRunFileGridDeformation.F90
	$(FC) $(CFLAGS) Runfiles/MainRunFileGridDeformation.F90 

## Goat.o 			: main run file goat
$(BUILDDIR)/goat.o: Runfiles/Goat.F90
	$(FC) $(CFLAGS) Runfiles/Goat.F90 -I$(BUILDDIR)

## Tests.o 			: all tests
$(BUILDDIR)/tests.o: Runfiles/Tests.F90 
	$(FC) $(CFLAGS) Runfiles/Tests.F90 -I$(BUILDDIR)

## Testc.o 			: C layer tests
$(BUILDDIR)/testc.o: Runfiles/Testc.c 
	$(CC) $(CCFLAGS) Runfiles/Testc.c -lcxsparse -I $(SUITESPARSEPATH) -I src/Clayer/Include -I$(BUILDDIR)

## goattranslator.o	: translator
$(BUILDDIR)/goattranslator.o: Runfiles/TranslateGOAToptionsFile.F90
	$(FC) $(CFLAGS) Runfiles/TranslateGOAToptionsFile.F90  -I$(BUILDDIR)

## shapeopt.o		: shape optimization 
$(BUILDDIR)/shapeopt.o : Runfiles/ShapeOptimization.F90 
	$(FC) $(CFLAGS) Runfiles/ShapeOptimization.F90 -I$(BUILDDIR)
	$(FC) $(CFLAGS) Runfiles/ShapeOptimization.F90 -I$(BUILDDIR)

## shapeopt_solps.o		: shape optimization with solps
$(BUILDDIR)/shapeopt_solps.o : Runfiles/ShapeOptimization.F90 
	$(FC) $(CFLAGS) Runfiles/ShapeOptimization.F90 -I$(BUILDDIR)

##
## % Folder compilation targets
## %===========================
## Constants 		: compile files containing constants
$(BUILDDIR)/Constants: $(CONSTANTS_FILES) 
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Constants 

## General 			: compile general files and modules
$(BUILDDIR)/General: $(GENERAL_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/General
	
## Modules			: compile modules 
$(BUILDDIR)/Modules: $(MODULE_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Modules

## Auxiliary			: compile auxiliary routines
$(BUILDDIR)/Auxiliary: $(AUXILIARY_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Auxiliary

## Drivers			: compile all (goat) driver routines
$(BUILDDIR)/Drivers: $(DRIVER_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Drivers

## SODrivers			: compile all shape optimization driver routines
$(BUILDDIR)/SODrivers: $(SODRIVER_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/SODrivers

## Setup			: compile all setup routines
$(BUILDDIR)/Setup: $(SETUP_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Setup

## IO_output 			: compile output routines
$(BUILDDIR)/IO_output: $(OUTPUT_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/IO_output

## IO_input			: compile input routines
$(BUILDDIR)/IO_input: $(INPUT_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/IO_input

## IO_b25			: compile b25 routines
$(BUILDDIR)/IO_b25: $(B25_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/IO_b25

## IO_carre			: compile carre routines
$(BUILDDIR)/IO_carre: $(CARRE_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/IO_carre

## Numerics		 	: compile numerics routines
$(BUILDDIR)/Numerics: $(NUMERICS_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Numerics

## Optimization		 	: compile optimization routines
$(BUILDDIR)/Optimization: $(OPTIMIZATION_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/Optimization

## Clayer 				: compile C interlayer routines
$(BUILDDIR)/Clayer: $(CLAYER_FILES)
	$(CC) $(CCFLAGS) $^ -I $(SUITESPARSEPATH) -I src/Clayer/Include -I$(BUILDDIR)
	touch $(BUILDDIR)/Clayer

## ClayerF 				: compile fortran modules in Clayer
$(BUILDDIR)/ClayerF: $(CLAYERF_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/ClayerF

## ShapeOptimization 			: compile shape optimization modules
$(BUILDDIR)/ShapeOptimization: $(SHAPEOPT_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/ShapeOptimization

## ShapeOptimizationSolps 			: compile shape optimization modules for SOLPS
$(BUILDDIR)/ShapeOptimizationSolps: $(SHAPEOPTSOLPS_FILES)
	$(FC) $(CFLAGS) -I$(B25LIBPATH) -I$(B25LIBPATH) -L$(B25LIBPATH) -l:libb2.a -lnetcdf $^ 
	touch $(BUILDDIR)/ShapeOptimizationSolps 


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

## cleanbuilds   		: clean builds directory (removes all subdirectories)
.PHONY: cleanbuilds 
cleanbuilds: 
	rm -r ./builds/*

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
