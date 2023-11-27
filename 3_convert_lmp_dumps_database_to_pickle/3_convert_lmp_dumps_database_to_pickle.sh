#!/bin/bash

function Convert_dumps_to_pickle {
    cs=$1
    specie=$2
    rcut=$3
    twojmax=$4
    nnn=$5
    delta=$6
    cstype=$7
    
    cp -r ../2_read_finite_temperature_database_apply_F_compute_descriptor/bnnn_${specie}_${cs}_rcut_${rcut}_nnn_${nnn}_delta_${delta}/dir.dumps/* bnnn_full_dumps_database

    python3 convert_dumps_to_pickle.py ${specie} ${cs} ${nnn} ${twojmax} ${cstype}
}

##########################################################################################################################################################################################################
# 3rd step :

rm -rf bnnn*

rcut=9
twojmax=8
nnn=24
delta=0.25

mkdir bnnn_full_dumps_database
mkdir bnnn_full_pickle_database

Convert_dumps_to_pickle bcc Fe ${rcut} ${twojmax} ${nnn} ${delta} 0
Convert_dumps_to_pickle fcc Al ${rcut} ${twojmax} ${nnn} ${delta} 1
Convert_dumps_to_pickle hcp Zr ${rcut} ${twojmax} ${nnn} ${delta} 2
Convert_dumps_to_pickle dia Si ${rcut} ${twojmax} ${nnn} ${delta} 3

##########################################################################################################################################################################################################
