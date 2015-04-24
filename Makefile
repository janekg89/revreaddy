CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
BUILDDIR := build
TARGET := bin/runner
 
SRCEXT := cpp
#SOURCES := src/main.cpp src/Random.cpp src/Particle.cpp  src/Simulation.cpp  src/Potential.cpp src/Trajectory.cpp src/Observable.cpp
SOURCES := $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
# objects without main.o
OBJECTSWM := $(subst build/main.o,,$(OBJECTS))
CFLAGS := -g -std=c++11#-Wall
LIB := -lgsl -lgslcblas -lm -lhdf5 -lhdf5_hl#-pthread -lmongoclient -L lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt
INC := -I include
UNITTESTS := $(shell find tests/unittests -type f -name *.h)
UNITTESTTARGET := bin/unittest
UNITTESTDIR := tests/unittests

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo "$(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET) $(UNITTESTTARGET)"; $(RM) -r $(BUILDDIR) $(TARGET) $(UNITTESTTARGET)

unittest: $(UNITTESTTARGET)

# Link test
$(UNITTESTTARGET): $(BUILDDIR)/unittest.o $(OBJECTSWM)
	@echo " Linking unittests...";
	@echo " $(CC) $^ -o $(UNITTESTTARGET) $(LIB)"; $(CC) $^ -o $(UNITTESTTARGET) $(LIB)

# How to build the unittest runner
$(BUILDDIR)/unittest.o: $(UNITTESTDIR)/unittest.cpp 
	@mkdir -p $(BUILDDIR)
	@echo " Build unittest...";
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $< "; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

# How to generate the test
$(UNITTESTDIR)/unittest.cpp: $(UNITTESTS)
	cxxtestgen --error-printer -o $@ $^

.PHONY: clean
