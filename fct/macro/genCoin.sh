#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage' ./genCoin.sh [emlist]"
    echo
}

plotany=0
single=0
while :; do
    if [ ${1}_ = "-h_" ]; then
        show_usage
        exit
    elif [ ${1}_ = "-s_" ]; then
        single=1
        shift 1
    elif [ ${1}_ = "-a_" ]; then
        plotany=1
        shift 1
    else
        break
    fi
done
list=${1}

if [ -z ${list} ]; then
    show_usage
    exit
elif [ ! -f "${list}" ]; then
    echo "emlist is not found"
    exit
fi

emcount=-1
boards=()
filenames=()

while read line; do
    echo "${line}"
    if [ -z "${line}" ]; then
        continue
    fi

    array=($(echo ${line}))
    this_emcount=${array[0]}
    this_board=${array[1]}
    this_date=${array[2]}
    this_filename=${array[3]}
    # echo "this_emcount  = ${this_emcount}"
    # echo "this_board    = ${this_board}"
    # echo "this_date     = ${this_date}"
    # echo "this_filename = ${this_filename}"

    this_root_dirname=$(dirname ${this_filename/data/ana})
    this_root_filename=$(basename ${this_filename/.dat/.root})
    # echo "this_root_dirname  = ${this_root_dirname}"
    # echo "this_root_filename = ${this_root_filename}"

    if [ ${emcount} -ne ${this_emcount} ]; then
        if [ ${#filenames[@]} -eq 8 ] || [ ${plotany} -eq 1 ] && [ ${#filenames[@]} -ne 0 ]; then
            marged_dirname=$(dirname ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
            marged_filename=$(basename ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
            marged_filename=${marged_filename/.root/.pdf}
            # echo "marged_dirname  = ${marged_dirname}"
            # echo "marged_filename = ${marged_filename}"

            if [ ! -d ${marged_dirname} ]; then
                mkdir -p ${marged_dirname}
            fi

            ${SOURCEDIR}/../build/genCoin \
                        ${SOURCEDIR}/../conf/genHist.conf \
                        $(echo ${boards[@]} | tr " " ",") \
                        $(echo ${filenames[@]} | tr " " ",") \
                        -o ${marged_dirname}/${marged_filename}

            if [ ${single} -eq 1 ]; then
                exit
            fi
        fi
        boards=()
        filenames=()
    fi

    emcount=${this_emcount}
    boards+=(${this_board})
    filenames+=(${this_root_dirname}/${this_root_filename})
done < <(cat ${list} | sort -k 1n,1 -k 2n,2)

if [ ${#filenames[@]} -eq 8 ] || [ ${plotany} -eq 1 ] && [ ${#filenames[@]} -ne 0 ]; then
    marged_dirname=$(dirname ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
    marged_filename=$(basename ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
    marged_filename=${marged_filename/.root/.pdf}
    # echo "marged_dirname  = ${marged_dirname}"
    # echo "marged_filename = ${marged_filename}"

    if [ ! -d ${marged_dirname} ]; then
        mkdir -p ${marged_dirname}
    fi

    ${SOURCEDIR}/../build/genCoin \
                ${SOURCEDIR}/../conf/genHist.conf \
                $(echo ${boards[@]} | tr " " ",") \
                $(echo ${filenames[@]} | tr " " ",") \
                -o ${marged_dirname}/${marged_filename}
fi
