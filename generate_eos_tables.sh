#!/bin/sh
# Generate all EOS tables for ANEOSmaterial.

MANEOS_EXE="../ANEOSmaterial/writeMANEOStable"
DATE=$(date +%d.%m.%Y)

#make clean

if [ ! -e "$MANEOS_EXE" ]; then
    # Generate M-ANEOS tables
    echo "Binary $MANEOS_EXE does not exist."
    if make "$MANEOS_EXE" 1> /dev/null 2>&1; then
        echo "Compiled writeMANEOStable."
    else
        echo "Compiling writeMANEOStable failed."
        exit 1
    fi
fi

# M-ANEOS iron alloy (Stewart 2020, http://doi.org/10.5281/zenodo.3866550)
INPUT="../maneos_input_files/Fe85Si15.input"
RHO0="7.51"
OUTPUT="MANEOStable_iron_Fe85Si15.in"
MAT_STR="M-ANEOS iron alloy Fe85Si15 ($DATE, Stewart 2020: http://doi.org/10.5281/zenodo.3866550)"

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS iron (Stewart 2020, http://doi.org/10.5281/zenodo.3866507)
INPUT="../maneos_input_files/iron.input"
RHO0="8.06"
OUTPUT="MANEOStable_iron.in"
MAT_STR="M-ANEOS iron ($DATE, Stewart 2020: http://doi.org/10.5281/zenodo.3866507)"

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS forsterite (Stewart 2019, http://doi.org/10.5281/zenodo.3478631)
INPUT="../maneos_input_files/forsterite.input"
RHO0="3.22"
OUTPUT="MANEOStable_forsterite.in"
MAT_STR="M-ANEOS forsterite ($DATE, Stewart 2019: http://doi.org/10.5281/zenodo.3478631)"

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS quartz (Stewart 2019, http://doi.wiley.com/10.1111/j.1945-5100.2007.tb01009.x)
INPUT="../maneos_input_files/quartz.input"
RHO0="2.65"
OUTPUT="MANEOStable_quartz.in"
MAT_STR="M-ANEOS quartz ($DATE, Melosh 2007: http://doi.wiley.com/10.1111/j.1945-5100.2007.tb01009.x)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS serpentine (Brookshaw 1998)
INPUT="../maneos_input_files/serpentine.input"
RHO0="2.50"
OUTPUT="MANEOStable_serpentine.in"
MAT_STR="M-ANEOS serpentine ($DATE, Brookshaw 1998)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS water (Stewart 2024, https://zenodo.org/records/14226694)
INPUT="../maneos_input_files/water.input"
RHO0="1.25"
OUTPUT="MANEOStable_water.in"
MAT_STR="M-ANEOS water ($DATE, Melosh 2007: https://zenodo.org/records/14226694)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi