objects = aneos.o 
mylibobjects = ANEOSmaterial.o interpBilinear.o

fortran_objects = libaneos.o

execs = testANEOSmaterial writeANEOStable writePhase

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas
FC := gfortran

CFLAGS ?= -O3 -Wall -std=c99

FFLAGS ?= $(CFLAGS)
LIBS ?= -lm -lgfortran $(GSL_LIB)

default: 
	@echo "Please specify which tool you want to make."	

all: default

testANEOSmaterial: testANEOSmaterial.o $(mylibobjects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	
writeANEOStable: writeANEOStable.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

writePhase: writePhase.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(execs) $(objects) $(fortran_objects) $(mylibobjects)

cleanall:
	rm -f $(execs) aneos.output fort.22 *.txt *.mat *.o *.pdf
