#!/bin/bash

hv=${1:-0}

if [ "${hv}" = "on" ]; then
    set -x

    ./setHV.py 1 57
    ./setHV.py 2 57
    ./setHV.py 3 57
    ./setHV.py 4 57
    ./setHV.py 5 57
    ./setHV.py 6 57
    ./setHV.py 7 57
    ./setHV.py 8 56

    set +x
elif [ "${hv}" = "off" ]; then
    set -x

    ./setHV.py 1 0
    ./setHV.py 2 0
    ./setHV.py 3 0
    ./setHV.py 4 0
    ./setHV.py 5 0
    ./setHV.py 6 0
    ./setHV.py 7 0
    ./setHV.py 8 0

    set +x
else
    set -x

    ./setHV.py 1 ${hv}
    ./setHV.py 2 ${hv}
    ./setHV.py 3 ${hv}
    ./setHV.py 4 ${hv}
    ./setHV.py 5 ${hv}
    ./setHV.py 6 ${hv}
    ./setHV.py 7 ${hv}
    ./setHV.py 8 ${hv}

    set +x
fi
