objects = aneos.o 
mylibobjects = ANEOSmaterial.o

defs = -DTILL_OUTPUT_ALL_WARNINGS -DTILL_VERBOSE

fortran_objects = libaneos.o

execs = testANEOSmaterial writeANEOStable

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas
FC := gfortran

CFLAGS ?= -O3

FFLAGS ?= $(CFLAGS)
LIBS ?= -lm -lgfortran $(GSL_LIB)

default: 
	@echo "Please specify which tool you want to make."	

all: default

testANEOSmaterial: testANEOSmaterial.o $(mylibobjects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	
writeANEOStable: writeANEOStable.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(execs) $(objects) $(fortran_objects) $(mylibobjects)

cleanall:
	rm -f $(execs) aneos.output fort.22 *.txt *.mat *.o *.pdf
