#!/bin/bash

vth=${1:-0}

if [ "${vth}" = "default" ]; then
    vth=100
fi

set -x

./setVth.py 1 $(echo "${vth} - 31.04" | bc)
./setVth.py 2 $(echo "${vth} - 70.54" | bc)
./setVth.py 3 $(echo "${vth} - 31.03" | bc)
./setVth.py 4 $(echo "${vth} - 31.04" | bc)
./setVth.py 5 $(echo "${vth} - 31.04" | bc)
./setVth.py 6 $(echo "${vth} - 19.09" | bc)
./setVth.py 7 $(echo "${vth} - 47.79" | bc)
./setVth.py 8 $(echo "${vth} - 27.63" | bc)

set +x
