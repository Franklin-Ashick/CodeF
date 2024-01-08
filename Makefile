# Compiler settings
GCC = gcc
GNU_MPI_CC = mpicc
INTEL_MPI_CC = mpiicc

# Compiler flags
CFLAGS = -fopenmp -std=c99
MPI_CFLAGS = -std=c99 # No OpenMP flag for MPI compilation

# Linker flags
LDFLAGS = -lm

# Source files
SOURCES = coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c
OMP_ONLY_SOURCE = main-openmp-only.c
MPI_SOURCE = main-mpi.c

# Default target
all: gomp-only gcomplete

# OpenMP Only Target
gomp-only: $(OMP_ONLY_SOURCE) $(SOURCES)
	$(GCC) $(CFLAGS) $(OMP_ONLY_SOURCE) $(SOURCES) -o gomp-only $(LDFLAGS)

# MPI Targets
gcomplete: $(MPI_SOURCE) $(SOURCES)
	$(GCC) $(CFLAGS) -I/opt/homebrew/include -L/opt/homebrew/lib -lmpi $(MPI_SOURCE) $(SOURCES) -o gcomplete $(LDFLAGS)

# Intel MPI Compiler Targets (if needed)
iserial: $(OMP_ONLY_SOURCE) $(SOURCES)
	$(INTEL_MPI_CC) $(CFLAGS) $(OMP_ONLY_SOURCE) $(SOURCES) -o iserial $(LDFLAGS)

icomplete: $(MPI_SOURCE) $(SOURCES)
	$(INTEL_MPI_CC) $(CFLAGS) $(MPI_SOURCE) $(SOURCES) -o icomplete $(LDFLAGS)

# Debug Target
debug: $(OMP_ONLY_SOURCE) $(SOURCES)
	$(GCC) -g $(CFLAGS) $(OMP_ONLY_SOURCE) $(SOURCES) -o debug $(LDFLAGS)

# Clean command
clean:
	rm -f *.exe *.o gomp-only gcomplete iserial icomplete debug

# Other custom targets can be added here
