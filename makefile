#Set the compilers and linker
FC=gfortran
CC=gcc
LD=gfortran

#Set the objects
# OBJS=			    	   precision_set.o \
# 									 read_gem_file.o \
# 									 maths_class.o \
# 									 error_class.o \
# 									 line_class.o \
# 									 material_class.o \
# 									 region_class.o \
# 									 region_class_1d.o \
# 									 timer_class.o \
# 									 matrix_class.o \
# 									 populate_class.o \
# 									 compressed_matrix_class.o \
# 									 tridiagonal_matrix_class.o \
# 									 populate_tridiagonal_class.o \
# 									 populate_compressed_class.o \
# 									 initialise_variables.o \
# 									 solver_class.o \
# 				    test_periodic.o
OBJS=			    	   precision_set.o \
									 maths_class.o \
									 line_class.o \
									 material_class.o \
									 region_class.o \
									 mesh_class.o \
									 region_class_1d.o \
									 region_class_2d.o \
									 matrix_class.o \
									 compressed_matrix_class.o \
									 tridiagonal_matrix_class.o \
									 populate_class.o \
									 read_gem_file.o \
									 error_class.o \
									 populate_tridiagonal_class.o \
									 populate_compressed_class.o \
									 populate_compressed_class_2d.o \
									 initialise_variables.o \
									 vtk_class.o \
									 solver_class.o \
									 timer_class.o \
				    test_2d.o


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fortran Compiler options
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FC_FLAGS    =           -fcheck=bounds -ffree-line-length-800 -O -pg -g -fimplicit-none -fbacktrace
#use O3 for optimising


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE MODS: OBJS .o REPLACED BY .mod
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MODS= $(OBJS:.o=.mod)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EXEC=test_2d

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
