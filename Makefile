# Executable
EXE = program.exe

# Git Repositary
GIT_URL = git@github.com:mazinmmustafa/Moments_Method_3D_Version1.git master

# Compilers
CXX = g++
CC = gcc
FC = gfortran

# Compiler Flags
CXXFLG = -Wall -Wextra -pedantic -std=c++14  
CFLG = -Wall -Wextra -pedantic -std=c17
F90FLG = -Wall -Wextra -pedantic
F77FLG =  

CLIB = -lgfortran -lquadmath -lm

# Optimization
OPT_LEVEL = -O2
OPT = $(OPT_LEVEL) -march=native -mtune=native
CXXOPT = $(OPT) 
COPT = $(OPT) 
F90OPT = $(OPT) 
F77OPT = -O2 
FLTO = -flto

# Directories
BDIR = bin
SDIR = src
ODIR = .obj
DDIR = .dep
HDIR = include

# Variables
CXXSRC = $(wildcard $(SDIR)/*.cpp)
CXXOBJ = $(patsubst $(SDIR)/%.cpp, $(ODIR)/%.o, $(CXXSRC))
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))
F90SRC = $(wildcard $(SDIR)/*.f90)
F90OBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(F90SRC))

CXXDEP = $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $(CXXOBJ))
CDEP = $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $(COBJ))

CXXDEPFLG = -MMD -MF $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $@) 

# Targets
.PHONY: build run clean clean_mesh clean_all figures

build: $(BDIR) $(ODIR) $(DDIR) $(BDIR)/$(EXE) 

run: build
	@echo executing $^:
	@./$(BDIR)/$(EXE)

$(BDIR): 
	@mkdir -pv $@

$(ODIR): 
	@mkdir -pv $@

$(DDIR): 
	@mkdir -pv $@

$(BDIR)/$(EXE): $(CXXOBJ) $(COBJ) $(F90OBJ) $(F77OBJ)
	$(CXX) $(FLTO) -o $@ $^ -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) -g $(FLTO) $(CXXFLG) $(CXXOPT) $(CXXDEPFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) -g $(FLTO) $(CFLG) $(COPT) $(CDEPFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -g $(FLTO) $(F90FLG) $(F90OPT) -o $@ -c $<

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(F77FLG) $(F77OPT) -o $@ -c $<

clean:
	@$(RM) -rv $(BDIR)/$(EXE) $(ODIR)/*.o $(DDIR)/*.d 

clean_mesh:
	@$(RM) -rv mesh/basis/* mesh/mesh/* mesh/shape.vtk

clean_all: clean clean_mesh
	@$(RM) -rv $(BDIR) $(ODIR) $(DDIR)

clean_data:
	@find ./data/ -type f -name '*.bin' | xargs $(RM) -rv
	@find ./data/ -type f -name '*.txt' | xargs $(RM) -rv
	@find ./figures/ -type f -name '*.pdf' | xargs $(RM) -rv
	@find ./figures/ -type f -name '*.txt' | xargs $(RM) -rv
	@find ./figures/ -type f -name '*.log' | xargs $(RM) -rv
	@find ./figures/ -type f -name '*.aux' | xargs $(RM) -rv
	@find ./figures/ -type f -name '*.tex' | xargs $(RM) -rv

.PHONY: valgrind 

valgrind: build
	valgrind --leak-check=full --track-origins=yes ./$(BDIR)/$(EXE) 

.PHONY: git_push git_pull gmsh

git_push: clean_all clean_data
	git add .
	git commit -m 'update' 
	git push --force --set-upstream $(GIT_URL)

git_pull: clean_all
	git pull $(GIT_URL)

gmsh:
	@python3 mesh/generate_mesh.py

figures:
	$(MAKE) -C figures/RCS1/

# Includes
-include $(CXXDEP) $(CDEP)

