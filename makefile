#Set the compilers and linker
FC=gfortran
CC=gcc
LD=gfortran

#Set the objects
OBJS=			    	   precision_set.o \
									 line_class.o \
									 material_class.o \
									 region_class.o \
									 region_class_1d.o \
									 timer_class.o \
									 read_gem_file.o \
									 error_class.o \
									 maths_class.o \
									 matrix_class.o \
									 nuclear_matrix_class.o \
									 initialise_variables.o \
				    test_eigenvalue.o
# OBJS=			    	   precision_set.o \
# 									 line_class.o \
# 									 material_class.o \
# 									 region_class.o \
# 									 region_class_1d.o \
# 									 timer_class.o \
# 									 read_gem_file.o \
# 									 error_class.o \
# 									 maths_class.o \
# 									 matrix_class.o \
# 									 nuclear_matrix_class.o \
# 									 initialise_variables.o \
# 				    test_regions.o
# OBJS=							 precision_set.o \
# 									 line_class.o \
# 									 material_class.o \
# 									 region_class.o \
# 									 region_class_1d.o \
# 									 timer_class.o \
# 									 read_gem_file.o \
# 									 error_class.o \
# 									 maths_class.o \
# 									 matrix_class.o \
# 									 nuclear_matrix_class.o \
# 									 initialise_variables.o \
# 				    test_classes.o


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fortran Compiler options
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FC_FLAGS    =           -fcheck=bounds -ffree-line-length-800 -O -pg
#use O3 for optimising


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE MODS: OBJS .o REPLACED BY .mod
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MODS= $(OBJS:.o=.mod)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EXEC=test_eigenvalue

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAKEFILE VARIABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEFAULT=makefile

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAKE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
all : $(EXEC)
debug : $(EXEC)
home: $(EXEC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LINK THE SOURCE (.o) INTO THE TARGET (EXEC) - explicit rule
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$(EXEC): $(OBJS)
	$(LD) $(FC_FLAGS)  -o $@ $^

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# COMPILE INFERRED SOURCE (.f90) INTO TARGET (.o) - inference rule
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%.o : %.f90
	$(FC) $(FC_FLAGS) -cpp -dM -c $< -o $@

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# UPDATE OBJS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%.f90 : $(DEFAULT)

#Clean the directory and any directories searched
.PHONY : clean
clean :
	rm -rf $(OBJS) $(EXEC) $(MODS) *.o *.mod
