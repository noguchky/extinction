#!/bin/bash

vth=${1:-0}

offset=(
    31.04
    70.54
    31.03
    31.04
    31.04
    19.09
    47.79
    27.63
)

if [ "${vth}" = "default" ]; then
    set -x

    ./setVth.py 1 $(echo "60.00 - ${offset[0]}" | bc)
    ./setVth.py 2 $(echo "60.00 - ${offset[1]}" | bc)
    ./setVth.py 3 $(echo "60.00 - ${offset[2]}" | bc)
    ./setVth.py 4 $(echo "60.00 - ${offset[3]}" | bc)
    ./setVth.py 5 $(echo "60.00 - ${offset[4]}" | bc)
    ./setVth.py 6 $(echo "60.00 - ${offset[5]}" | bc)
    ./setVth.py 7 $(echo "60.00 - ${offset[6]}" | bc)
    ./setVth.py 8 $(echo "60.00 - ${offset[7]}" | bc)

    set +x
else
    set -x

    ./setVth.py 1 $(echo "${vth} - ${offset[0]}" | bc)
    ./setVth.py 2 $(echo "${vth} - ${offset[1]}" | bc)
    ./setVth.py 3 $(echo "${vth} - ${offset[2]}" | bc)
    ./setVth.py 4 $(echo "${vth} - ${offset[3]}" | bc)
    ./setVth.py 5 $(echo "${vth} - ${offset[4]}" | bc)
    ./setVth.py 6 $(echo "${vth} - ${offset[5]}" | bc)
    ./setVth.py 7 $(echo "${vth} - ${offset[6]}" | bc)
    ./setVth.py 8 $(echo "${vth} - ${offset[7]}" | bc)

    set +x
fi
