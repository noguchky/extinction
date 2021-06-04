#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage' ./genHist.sh [emlist]"
    echo
}

confFilename="${SOURCEDIR}/../conf/genHist.conf"
single=0
plotany=0
delayopt=
while :; do
    if   [ ${1}_ = "-h_" ]; then
        show_usage
        exit
    elif [ ${1}_ = "-c_" ]; then
        confFilename=${2}
        shift 2
    elif [ ${1}_ = "-s_" ]; then
        single=1
        shift 1
    elif [ ${1}_ = "-a_" ]; then
        plotany=1
        shift 1
    elif [ ${1}_ = "-d_" ]; then
        delayopt="-d"
        shift 1
    else
        break
    fi
done
list=${@}

if [ -z "${list}" ]; then
    show_usage
    exit
fi

for fname in ${list}; do
    if [ ! -f "${fname}" ]; then
        echo "emlist\"${fname}\" is not found"
        exit
    fi
done

emcount=-1
boards=()
filenames=()

function exec_genHist() {
    if [ ${#filenames[@]} -eq 8 ] || [ ${plotany} -eq 1 ] && [ ${#filenames[@]} -ne 0 ]; then
        marged_dirname=$(dirname ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
        marged_filename=$(basename ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
        marged_filename=${marged_filename/.root/.pdf}
        # echo "marged_dirname  = ${marged_dirname}"
        # echo "marged_filename = ${marged_filename}"

        if [ ! -d ${marged_dirname} ]; then
            mkdir -p ${marged_dirname}
        fi

        ${SOURCEDIR}/../build/genHist \
                    ${confFilename} \
                    $(echo ${boards[@]} | tr " " ",") \
                    $(echo ${filenames[@]} | tr " " ",") \
                    -o ${marged_dirname}/${marged_filename} \
                    ${delayopt}

        if [ ${single} -eq 1 ]; then
            exit
        fi
    fi
}

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

    if [ ! -f ${this_root_dirname}/${this_root_filename} ]; then
        echo "root file is not exist, need to decode"
        exit
    fi

    if [ ${emcount} -ne ${this_emcount} ]; then
        exec_genHist
        boards=()
        filenames=()
    fi

    emcount=${this_emcount}
    boards+=(${this_board})
    filenames+=(${this_root_dirname}/${this_root_filename})
done < <(cat ${list} | sort -k 1n,1 -k 2n,2)

exec_genHist
