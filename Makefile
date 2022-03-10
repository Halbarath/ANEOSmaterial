objects = aneos.o 

mylibobjects = ANEOSmaterial.o interpBilinear.o

# Compile with ANEOS
#fortran_objects = libaneos.o
# Compile with M-ANEOS
aneos_path = ../M-ANEOS/src
libmaneos_objects = $(aneos_path)/ANDATA.o \
                  $(aneos_path)/ANDEBY.o \
                  $(aneos_path)/ANE2PH.o \
                  $(aneos_path)/ANEHPP.o \
                  $(aneos_path)/ANEI3.o \
                  $(aneos_path)/ANEINI.o \
                  $(aneos_path)/ANELSM.o \
                  $(aneos_path)/ANEOS1.o \
	              $(aneos_path)/ANEOS2.o \
	              $(aneos_path)/ANEOSD.o \
	              $(aneos_path)/ANEOS.o \
	              $(aneos_path)/ANEOSI.o \
	              $(aneos_path)/ANEOSS.o \
	              $(aneos_path)/ANEOSV.o \
	              $(aneos_path)/ANEVAL.o \
	              $(aneos_path)/ANHUG.o \
	              $(aneos_path)/ANION1.o \
	              $(aneos_path)/ANION2.o \
	              $(aneos_path)/ANMARK.o \
	              $(aneos_path)/ANMAXW.o \
	              $(aneos_path)/ANN1AS.o \
	              $(aneos_path)/ANN1VL.o \
	              $(aneos_path)/ANN2AS.o \
	              $(aneos_path)/ANNDPR.o \
	              $(aneos_path)/ANPHAS.o \
	              $(aneos_path)/ANPHTR.o \
	              $(aneos_path)/ANPRTR.o \
	              $(aneos_path)/ANSIOF.o \
	              $(aneos_path)/ANSMFT.o \
	              $(aneos_path)/ANTOMF.o \
	              $(aneos_path)/ANUEOS.o \
	              $(aneos_path)/ANUSET.o \
	              $(aneos_path)/ANWARN.o \
	              $(aneos_path)/ANZRTR.o
fortran_objects = $(libmaneos_objects) $(aneos_path)/ANEOSINIT.o

execs = testANEOSmaterial writeANEOStable writeMANEOStable writePhase writePressureTable aneoscall tipsy_iphase_array aneostableinfo tipsy_ascii

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas
FC := gfortran

CFLAGS ?= -O3 -Wall -std=c99

#FFLAGS ?= $(CFLAGS) -std=legacy
FFLAGS = -O3
LIBS ?= -lm -lgfortran $(GSL_LIB)

default: 
	@echo "Please specify which tool you want to make."	

all: default

testANEOSmaterial: testANEOSmaterial.o $(mylibobjects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testMANEOStable: testMANEOStable.o $(mylibobjects) $(fortran_objects) $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

testMANEOStableint: testMANEOStableint.o $(mylibobjects) $(fortran_objects) $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcPressureRhoT: calcPressureRhoT.o $(mylibobjects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

calcDensityPT: calcDensityPT.o $(mylibobjects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

writeANEOStable: writeANEOStable.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

writeMANEOStable: writeMANEOStable.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

writePhase: writePhase.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
    
writePhaseMANEOS: writePhaseMANEOS.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

writePressureTable: writePressureTable.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

aneoscall: aneoscall.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

# Direct call to M-ANEOS. Be sure to compile with the correct object files.
maneoscall: maneoscall.o $(objects) $(fortran_objects) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

tipsy_iphase_array: tipsy_iphase_array.o $(mylibobjects) tipsy.o 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

tipsy_ascii: tipsy_ascii.o $(objects) $(fortran_objects) tipsy.o 
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
    
aneostableinfo: aneostableinfo.o $(mylibobjects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(execs) $(objects) $(fortran_objects) $(mylibobjects) *.o

cleanall:
	rm -f $(execs) aneos.output fort.22 *.txt *.mat *.o *.pdf
