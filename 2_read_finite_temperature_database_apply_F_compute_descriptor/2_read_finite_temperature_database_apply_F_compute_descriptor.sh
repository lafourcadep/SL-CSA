#!/bin/bash

LMP_EXEC=/home/lafourcadep/CODES/LAMMPS/lammps_repo/build/lmp

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

##########################################################################################################################################################################################################
# 2nd step : read finite temperature restart files, apply deformation tensor and compute descriptor to build database in LAMMPS dump file format
rm -rf bnnn_*

rcut=9
twojmax=8
nnn=24
wmode=1
delta=0.25

Launch_BSO4_computation bcc Fe ${rcut} ${twojmax} ${nnn} ${wmode} ${delta}
Launch_BSO4_computation fcc Al ${rcut} ${twojmax} ${nnn} ${wmode} ${delta}
Launch_BSO4_computation hcp Zr ${rcut} ${twojmax} ${nnn} ${wmode} ${delta}
Launch_BSO4_computation dia Si ${rcut} ${twojmax} ${nnn} ${wmode} ${delta}

##########################################################################################################################################################################################################
