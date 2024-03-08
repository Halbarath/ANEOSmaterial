#!/bin/sh
# Generate all EOS tables for ANEOSmaterial.

MANEOS_EXE="writeMANEOStable"
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
INPUT="../M-ANEOS/input_eoslib/maneos_Fe85Si15_sts2020.input"
RHO0="7.51"
OUTPUT="MANEOStable_iron_Fe85Si15.in"
MAT_STR="M-ANEOS iron alloy Fe85Si15 ($DATE, Stewart 2020: http://doi.org/10.5281/zenodo.3866550)"

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS iron (Stewart 2020, http://doi.org/10.5281/zenodo.3866507)
INPUT="../M-ANEOS/input_eoslib/maneos_iron_sts2020.input"
RHO0="8.06"
OUTPUT="MANEOStable_iron.in"
MAT_STR="M-ANEOS iron ($DATE, Stewart 2020: http://doi.org/10.5281/zenodo.3866507)"

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS forsterite (Stewart 2019, http://doi.org/10.5281/zenodo.3478631)
INPUT="../M-ANEOS/input_eoslib/maneos_fosterite_sts2019.input"
RHO0="3.22"
OUTPUT="MANEOStable_fosterite.in"
MAT_STR="M-ANEOS forsterite ($DATE, Stewart 2019: http://doi.org/10.5281/zenodo.3478631)"

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi

# M-ANEOS forsterite (Stewart 2019, http://doi.org/10.5281/zenodo.3478631)
INPUT="../M-ANEOS/input_eoslib/maneos_dunite_gsc.input"
RHO0="2.65"
OUTPUT="MANEOStable_dunite.in"
MAT_STR="M-ANEOS dunite ($DATE, Collins and Melosh 2014)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1

if [ -e "$INPUT" ]; then
    "./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
else
    echo "ANEOS input file $INPUT not found."
fi
