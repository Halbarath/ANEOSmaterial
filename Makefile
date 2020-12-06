objects = aneos.o 

mylibobjects = ANEOSmaterial.o interpBilinear.o

# Compile with ANEOS
fortran_objects = libaneos.o
# Compile with MANEOS
#fortran_objects = libmaneos.o

execs = testANEOSmaterial writeANEOStable writePhase writePressureTable

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
    
writePressureTable: writePressureTable.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

aneoscall: aneoscall.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

aneoscalcentropy: aneoscalcentropy.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

tipsy_iphase_array: tipsy_iphase_array.o $(objects) $(fortran_objects) tipsy.o 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(execs) $(objects) $(fortran_objects) $(mylibobjects)

cleanall:
	rm -f $(execs) aneos.output fort.22 *.txt *.mat *.o *.pdf
