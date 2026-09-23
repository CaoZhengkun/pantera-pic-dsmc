############################################################
#                         MAKEFILE                         #
############################################################

# GNU Compiler
CMP  = mpifort -c -cpp -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LNK  = mpifort -cpp
OPTF = -O3 -fimplicit-none

SRCDIR = src/
BUILDDIR = src/

# A few more compiler flags for GNU fortran that could be useful:
#-O0 -Ofast -Wall -Wextra -Warray-temporaries -ggdb3 -pedantic -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Debug options
#-march=native -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Aggressive optimization options 

# Objects: list of all objects *.o
OBJS = $(BUILDDIR)mpi_common.o  $(BUILDDIR)velocity_distribution.o  $(BUILDDIR)global.o  $(BUILDDIR)screen.o  $(BUILDDIR)tools.o  $(BUILDDIR)field_dof_utils.o  $(BUILDDIR)external_field.o  $(BUILDDIR)initialization.o  $(BUILDDIR)timecycle.o  $(BUILDDIR)grid_and_partition.o  $(BUILDDIR)particle.o  $(BUILDDIR)collisions.o  $(BUILDDIR)postprocess.o  $(BUILDDIR)fields.o  $(BUILDDIR)mt19937.o  $(BUILDDIR)surface_emission_utils.o  $(BUILDDIR)periodic_mesh_utils.o  $(BUILDDIR)washboard.o

#  The following variable must either be a path to petsc.pc or just "petsc" if petsc.pc
#  has been installed to a system location or can be found in PKG_CONFIG_PATH.
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm
CUDAC := $(shell pkg-config --variable=cudacompiler $(PACKAGES))
CUDAC_FLAGS := $(shell pkg-config --variable=cudaflags_extra $(PACKAGES))
CUDA_LIB := $(shell pkg-config --variable=cudalib $(PACKAGES))
CUDA_INCLUDE := $(shell pkg-config --variable=cudainclude $(PACKAGES))


all: pantera.exe

debug: CMP += -g
debug: LNK += -g
debug: OPTF = -O0 -fimplicit-none -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
debug: pantera.exe

# Executable generation by the linker
pantera.exe: createbuilddir
pantera.exe: $(BUILDDIR)pantera.o $(OBJS) 
	$(LNK) $(OPTF) $(BUILDDIR)pantera.o $(OBJS) \
	            -o pantera.exe -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

stage0_utilities_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage0_utilities.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage0_utilities.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

.PHONY: test-stage0-utilities
test-stage0-utilities: stage0_utilities_test.exe
	./stage0_utilities_test.exe

stage1_source_region_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage1_source_region.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage1_source_region.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

.PHONY: test-stage1-source-region
test-stage1-source-region: stage1_source_region_test.exe
	./stage1_source_region_test.exe

stage3_periodic_mesh_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage3_periodic_mesh.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage3_periodic_mesh.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

.PHONY: test-stage3-periodic-mesh
test-stage3-periodic-mesh: stage3_periodic_mesh_test.exe
	./stage3_periodic_mesh_test.exe

stage3_periodic_particle_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage3_periodic_particle.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage3_periodic_particle.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

.PHONY: test-stage3-periodic-poisson
stage3_periodic_poisson_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage3_periodic_poisson.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage3_periodic_poisson.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

.PHONY: test-stage3-periodic-particle
test-stage3-periodic-particle: stage3_periodic_particle_test.exe
	./stage3_periodic_particle_test.exe

test-stage3-periodic-poisson: stage3_periodic_poisson_test.exe
	./stage3_periodic_poisson_test.exe

.PHONY: test-stage3-solver-extensions
stage3_solver_extensions_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage3_solver_extensions.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage3_solver_extensions.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

stage3_external_field_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage3_external_field.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage3_external_field.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

stage3_source_diagnostics_test.exe: createbuilddir $(OBJS) Verification/tests/test_stage3_source_diagnostics.f90
	$(LNK) $(OPTF) -I$(BUILDDIR) Verification/tests/test_stage3_source_diagnostics.f90 $(OBJS) \
	            -o $@ -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

test-stage3-solver-extensions: stage3_solver_extensions_test.exe stage3_external_field_test.exe
	./stage3_solver_extensions_test.exe
	./stage3_external_field_test.exe

test-stage3-source-diagnostics: stage3_source_diagnostics_test.exe
	./stage3_source_diagnostics_test.exe

# Objects generation
$(BUILDDIR)pantera.o: $(SRCDIR)pantera.f90  $(OBJS) createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)pantera.f90

$(BUILDDIR)global.o: $(SRCDIR)global.f90  $(BUILDDIR)mpi_common.o  $(BUILDDIR)particle.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)global.f90

$(BUILDDIR)timecycle.o: $(SRCDIR)timecycle.f90  $(BUILDDIR)global.o  $(BUILDDIR)particle.o  $(BUILDDIR)screen.o  $(BUILDDIR)velocity_distribution.o  $(BUILDDIR)collisions.o  $(BUILDDIR)postprocess.o  $(BUILDDIR)fields.o  $(BUILDDIR)washboard.o  $(BUILDDIR)grid_and_partition.o  $(BUILDDIR)surface_emission_utils.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)timecycle.f90

$(BUILDDIR)initialization.o: $(SRCDIR)initialization.f90  $(BUILDDIR)global.o  $(BUILDDIR)velocity_distribution.o  $(BUILDDIR)tools.o  $(BUILDDIR)surface_emission_utils.o  $(BUILDDIR)external_field.o  $(BUILDDIR)grid_and_partition.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)initialization.f90

$(BUILDDIR)tools.o: $(SRCDIR)tools.f90  $(BUILDDIR)mpi_common.o  $(BUILDDIR)global.o  $(BUILDDIR)screen.o  $(BUILDDIR)mt19937.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)tools.f90

$(BUILDDIR)field_dof_utils.o: $(SRCDIR)field_dof_utils.f90  $(BUILDDIR)global.o  $(BUILDDIR)screen.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)field_dof_utils.f90

$(BUILDDIR)external_field.o: $(SRCDIR)external_field.f90 createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)external_field.f90

$(BUILDDIR)grid_and_partition.o: $(SRCDIR)grid_and_partition.f90  $(BUILDDIR)mpi_common.o  $(BUILDDIR)global.o  $(BUILDDIR)tools.o  $(BUILDDIR)screen.o  $(BUILDDIR)periodic_mesh_utils.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)grid_and_partition.f90

$(BUILDDIR)screen.o: $(SRCDIR)screen.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)screen.f90

$(BUILDDIR)mpi_common.o: $(SRCDIR)mpi_common.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)mpi_common.f90

$(BUILDDIR)particle.o: $(SRCDIR)particle.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)particle.f90

$(BUILDDIR)velocity_distribution.o: $(SRCDIR)velocity_distribution.f90  $(BUILDDIR)tools.o  $(BUILDDIR)global.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)velocity_distribution.f90

$(BUILDDIR)collisions.o: $(SRCDIR)collisions.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)collisions.f90
	
$(BUILDDIR)postprocess.o: $(SRCDIR)postprocess.f90  $(BUILDDIR)fields.o  $(BUILDDIR)surface_emission_utils.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)postprocess.f90

$(BUILDDIR)fields.o: $(SRCDIR)fields.f90  $(BUILDDIR)global.o  $(BUILDDIR)screen.o  $(BUILDDIR)tools.o  $(BUILDDIR)field_dof_utils.o  $(BUILDDIR)external_field.o  $(BUILDDIR)grid_and_partition.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)fields.f90
	
$(BUILDDIR)washboard.o: $(SRCDIR)washboard.f90  $(BUILDDIR)tools.o  $(BUILDDIR)global.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)washboard.f90

$(BUILDDIR)mt19937.o: $(SRCDIR)mt19937.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)mt19937.f90

$(BUILDDIR)surface_emission_utils.o: $(SRCDIR)surface_emission_utils.f90  $(BUILDDIR)mt19937.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)surface_emission_utils.f90

$(BUILDDIR)periodic_mesh_utils.o: $(SRCDIR)periodic_mesh_utils.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)periodic_mesh_utils.f90

	
# Cleaning command
createbuilddir:  clean
	mkdir  build

clean: 
	@echo cleaning objects, modules and executables 
	rm -rf build *.exe src/*.mod src/*.o

cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*

