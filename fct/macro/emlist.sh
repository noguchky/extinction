#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage' ./genHist.sh [year] [month] (day) (hour) (month) (second)"
    echo
    echo "Options'" 
    echo " -e [chanel]        Set event-match channel"
    echo
}

EM_CH=27
while :; do
    if   [ ${1}_ = "-e_" ]; then
        EM_CH=${2}
        shift 2
    elif [ ${1}_ = "-h_" ]; then
        show_usage
        exit
    else
        break
    fi
done

year=${1}
month=${2}
day=${3:-xx}
hour=${4:-xx}
minute=${5:-xx}
second=${6:-xx}
tag=${7:-${year}${month}${day}${hour}${minute}${second}}

listname=./list/em_${tag}.txt

if [ -z ${year}${month} ]; then
    show_usage
    exit
fi

for board in {0..7}; do
    for filename in ../data/id00$(expr ${board} + 56)/${year}${month}${day//x/}/${hour//x/}*/fct_*_${year}${month}${day//x/}${hour//x/}${minute//x/}${second//x/}*.dat; do
        fullpath_filename=$(readlink -f ${filename})
        if [ -n "${fullpath_filename}" ]; then
            ${SOURCEDIR}/../build/emcount ${fullpath_filename} -b ${board} -e ${EM_CH}
        else
            echo "'${filename}' has no link" > /dev/stderr
        fi
    done
done | tee ${listname}

echo
echo "Info in emlist.sh: Event match list added to text file ${listname}"
echo
