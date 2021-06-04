#!/bin/bash

cd ~/noguchi/extinction/fct/macro

source ${SHOME}/common.sh

AHOME=${DHOME}/../ana/


unset dfiles
for ip in ${ips[@]}; do
    dir=$DHOME/$(printf "id%04d" $ip)/$(date '+%Y%m%d')
    dfname=$(find $dir -type f -name "$(printf "fct_id%04d" $ip)*" | sort | tail -n 1)
    echo $dfname
    if [ -z "$dfiles" ]; then
        dfiles=($dfname)
    else
        dfiles+=($dfname)
    fi
done

echo "${dfiles[@]}"

unset rfiles
unset odir
unset ofname
for dfname in ${dfiles[@]}; do
    rdir=$(dirname ${dfname/data/ana})
    rfname=${rdir}/$(basename ${dfname/.dat/.root})

    if [ ! -d ${rdir} ]; then
        mkdir -p ${rdir}
    fi

    if [ ! -f ${rfname} ]; then
        ../build/decoder ${dfname} -o ${rfname} &
    fi
    if [ -z "$rfiles" ]; then
        rfiles=($rfname)

        odir=$(dirname ${rfname//id[0-9][0-9][0-9][0-9]/marged})
        ofname=$(basename ${rfname//id[0-9][0-9][0-9][0-9]/marged})
        ofname=${odir}/${ofname/.root/.pdf}
    else
        rfiles+=($rfname)
    fi
done
wait

echo "${rfiles[@]}"
echo "${ofname}"

if [ ! -d ${odir} ]; then
    mkdir -p ${odir}
fi

../build/genHist \
    ../conf/genHist.conf \
    0,1,2,3,4,5,6,7 \
    $(echo ${rfiles[@]} | tr " " ",") \
    -o ${ofname}
