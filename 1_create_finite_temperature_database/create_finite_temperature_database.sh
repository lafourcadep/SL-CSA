#!/bin/bash

LMP_EXEC=/home/lafourcadep/CODES/LAMMPS/lammps_repo/build/lmp

function Launch_NPT_limited {
    cs=$1
    specie=$2
    temp=$3
    hrate=$4
    fractmelt=$5
    
    cp -r template npt_${specie}_${cs}

    cd npt_${specie}_${cs}/

    cp ../../0_heat_to_melt/npt_${specie}_${cs}/dir.data/low_T_low_P_equilibrated_${specie}_${cs}.data dir.structure/

    echo ${temp}
    sed -i "s/TEMPERATURE/${temp}/g"    dir.inputs/in_lammps_npt
    sed -i "s/FRACTMELT/${fractmelt}/g" dir.inputs/in_lammps_npt    
    sed -i "s/RATE/${hrate}/g"          dir.inputs/in_lammps_npt
    sed -i "s/SPECIE/${specie}/g"       dir.inputs/in_lammps_npt
    sed -i "s/STRUCTURE/${cs}/g"        dir.inputs/in_lammps_npt

    mpirun -np 8 ${LMP_EXEC} -in dir.inputs/in_lammps_npt
    
    cd ../
}

##########################################################################################################################################################################################################
# 1st step : run NPT temperature up to a fraction of melting temperature computed at 1st step

rm -rf npt_*

tmelt_Fe=2100.
tmelt_Al=1168.
tmelt_Zr=2600.
tmelt_Si=2597.
heatrate=20
fractmelt=67
Launch_NPT_limited bcc Fe ${tmelt_Fe} ${heatrate} ${fractmelt}
Launch_NPT_limited fcc Al ${tmelt_Al} ${heatrate} ${fractmelt}
Launch_NPT_limited hcp Zr ${tmelt_Zr} ${heatrate} ${fractmelt}
Launch_NPT_limited dia Si ${tmelt_Si} ${heatrate} ${fractmelt}

##########################################################################################################################################################################################################
