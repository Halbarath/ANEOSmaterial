#!/bin/sh
# Generate all EOS tables for ANEOSmaterial.

MANEOS_EXE="writeMANEOStable"
DATE=$(date +%d-%m-%Y)

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
INPUT="aneos-Fe85Si15-2020.input"
RHO0="7.51"
OUTPUT="MANEOStable_iron_Fe85Si15.in"
MAT_STR="M-ANEOS iron alloy Fe85Si15 ($DATE, Stewart 2020: http://doi.org/10.5281/zenodo.3866550)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1

# M-ANEOS iron (Stewart 2020, http://doi.org/10.5281/zenodo.3866507)
INPUT="aneos-iron-2020.input"
RHO0="8.06"
OUTPUT="MANEOStable_iron.in"
MAT_STR="M-ANEOS iron ($DATE, Stewart 2020: http://doi.org/10.5281/zenodo.3866507)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1

# M-ANEOS forsterite (Stewart 2019, http://doi.org/10.5281/zenodo.3478631)
INPUT="aneos-forsterite-2019.input"
RHO0="3.22"
OUTPUT="MANEOStable_forsterite.in"
MAT_STR="M-ANEOS forsterite ($DATE, Stewart 2019: http://doi.org/10.5281/zenodo.3478631)"

"./$MANEOS_EXE" "$INPUT" "$RHO0" "$OUTPUT" "$MAT_STR" 1> /dev/null 2>&1
