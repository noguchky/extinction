#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage: ./decode.sh (options) [emlist]"
    echo
    echo "Options:" 
    echo " -e [chanel]        Set event-match channel"
    echo
}

EM_CH=27
while :; do
    if   [ ${1}_ == "-e_" ]; then
        EM_CH=${2}
        shift 2
    elif [ ${1}_ == "-h_" ]; then
        show_usage
        exit
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

while read line; do
    echo ${line}
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
        wait
    fi

    if [ ! -d ${this_root_dirname} ]; then
        mkdir -p ${this_root_dirname}
    fi

    ${SOURCEDIR}/../build/decoder \
                ${this_filename} \
                -e ${EM_CH} \
                -o ${this_root_dirname}/${this_root_filename} &

    emcount=${this_emcount}
done < <(cat ${list} | sort -k 1n,1 -k 2n,2)
wait
