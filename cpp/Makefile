# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -O2

# PETSc settings
PETSC_DIR = /usr/lib/petscdir/petsc3.15
PETSC_ARCH = x86_64-linux-gnu-real

# PETSc library and include directories
PETSC_LIB = $(PETSC_DIR)/$(PETSC_ARCH)/lib
PETSC_INC = $(PETSC_DIR)/$(PETSC_ARCH)/include

# MPI include directory
MPI_INC = /usr/lib/x86_64-linux-gnu/openmpi/include

# Executable name
EXECUTABLE = your_program

# Source files
SRCS = Bubble_Growth_Modelling.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Compile the program
$(EXECUTABLE): $(OBJS)
	$(CXX) -o $@ $(OBJS) -L$(PETSC_LIB) -lpetsc

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(PETSC_INC) -I$(MPI_INC) -c $< -o $@

# Clean up
clean:
	rm -f $(EXECUTABLE) $(OBJS)

