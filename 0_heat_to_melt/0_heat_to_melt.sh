#!/bin/bash

LMP_EXEC=/home/lafourcadep/CODES/LAMMPS/lammps_repo/build/lmp

function Launch_NPT {
    cs=$1
    specie=$2
    temp=$3
    hrate=$4

    cp -r template npt_${specie}_${cs}

    cd npt_${specie}_${cs}/

    sed -i "s/TEMPERATURE/${temp}/g" dir.inputs/in_lammps_npt
    sed -i "s/RATE/${hrate}/g"       dir.inputs/in_lammps_npt
    sed -i "s/SPECIE/${specie}/g"    dir.inputs/in_lammps_npt
    sed -i "s/STRUCTURE/${cs}/g"     dir.inputs/in_lammps_npt

    mpirun -np 8 ${LMP_EXEC} -in dir.inputs/in_lammps_npt
    
    cd ../
}

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

function Launch_BSO4_computation {
    cs=$1
    specie=$2
    rcut=$3
    twojmax=$4
    nnn=$5
    wmode=$6
    delta=$7
    
    cp -r template bnnn_${specie}_${cs}_rcut_${rcut}_nnn_${nnn}_delta_${delta}
    
    cp -r ../1_create_finite_temperature_database/npt_${specie}_${cs}/dir.restarts bnnn_${specie}_${cs}_rcut_${rcut}_nnn_${nnn}_delta_${delta}/

    cd bnnn_${specie}_${cs}_rcut_${rcut}_nnn_${nnn}_delta_${delta}/

    sed -i "s/RCUT/${rcut}/g"         dir.inputs/${specie}.potential
    
    python3 generate_deformation_files.py

    cnt=0
    for entry in `ls dir.restarts/`;
    do
	timestep=$(echo $entry | awk -F[..] '{print $2}')
	echo ${timestep}
	cp dir.inputs/in_lammps_bnnn dir.inputs/in_lammps_bnnn_${timestep}

    	sed -i "s/TIMESTEP/${timestep}/g" dir.inputs/in_lammps_bnnn_${timestep}
    	sed -i "s/COUNTER/${cnt}/g"       dir.inputs/in_lammps_bnnn_${timestep}
    	sed -i "s/SPECIE/${specie}/g"     dir.inputs/in_lammps_bnnn_${timestep}
    	sed -i "s/STRUCTURE/${cs}/g"      dir.inputs/in_lammps_bnnn_${timestep}

	sed -i "s/RCUT/${rcut}/g"         dir.inputs/in_lammps_bnnn_${timestep}
	sed -i "s/TWOJMAX/${twojmax}/g"   dir.inputs/in_lammps_bnnn_${timestep}
    	sed -i "s/NNEIGH/${nnn}/g"        dir.inputs/in_lammps_bnnn_${timestep}
    	sed -i "s/WMODEVAL/${wmode}/g"    dir.inputs/in_lammps_bnnn_${timestep}
    	sed -i "s/DELTAVAL/${delta}/g"    dir.inputs/in_lammps_bnnn_${timestep}
	
	mpirun -np 8 ${LMP_EXEC} -in dir.inputs/in_lammps_bnnn_${timestep}
	
	let cnt=cnt+1
    done

    cd ../
    
}

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
# 0th step : run NPT temperature ramp up to tempmax to force material to melt (ideally tempmax should be large enough to reach melting)

rm -rf npt_*

# target temperature in K
tempmax=3000
# heating rate in K/ps
heatrate=20
Launch_NPT bcc Fe ${tempmax} ${heatrate}
Launch_NPT fcc Al ${tempmax} ${heatrate} 
Launch_NPT hcp Zr ${tempmax} ${heatrate}
Launch_NPT dia Si ${tempmax} ${heatrate}

# retrieve list of melting temperatures
result=`python3 identify_melting_temperature.py`
IFS=' '
read -a melting_ts <<< $result
tmelt_Fe=${melting_ts[0]}
tmelt_Al=${melting_ts[1]}
tmelt_Zr=${melting_ts[2]}
tmelt_Si=${melting_ts[3]}
