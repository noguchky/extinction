#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage' ./emlist.sh [year] [month] (day) (hour) (month) (second)"
    echo
    echo "Options'" 
    echo " -e [chanel]        Set event-match channel"
    echo
}

prefix=fct
boards=8
idoffset=56
emchannel=27
while :; do
    if   [ ${1}_ = "-e_" ]; then
        emchannel=${2}
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

if [ -z "${year}${month}" ]; then
    show_usage
    exit
fi

for ((board=0; board < boards; ++board)); do
    idstr=id00$(printf %02d $(expr ${board} + ${idoffset}))
    for filename in ../data/${idstr}/${year}${month}${day//x/}*/${hour//x/}*/${prefix}_${idstr}_${year}${month}${day//x/}${hour//x/}${minute//x/}${second//x/}*.dat; do
        fullpath_filename=$(echo ${SOURCEDIR}/${filename} | sed -r -e 's#/[^/]+/+\.\./#/#g' | sed -r -e 's#/+#/#g')
        if [ -n "${fullpath_filename}" ]; then
            ${SOURCEDIR}/../build/emcount ${fullpath_filename} -b ${board} -e ${emchannel}
        else
            echo "'${filename}' has no link" > /dev/stderr
        fi
    done
done | tee ${listname}

echo
echo "Info in emlist.sh: Event match list added to text file ${listname}"
echo
