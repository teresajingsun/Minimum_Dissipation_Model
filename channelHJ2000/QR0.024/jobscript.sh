#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30GB
#SBATCH --time=06-00:00:00
#SBATCH --job-name=QR2000

if [ -z "constant" ]
then 
    echo "ERROR: No input specified"
    exit -1
fi

if [ -e "constant" ]  
then  
    module purge
    module load OpenFOAM/v2006-foss-2020a
    source /software/software/OpenFOAM/v2006-foss-2020a/OpenFOAM-v2006/etc/bashrc

    echo "Processing LES: " "QR2000"

    ./Allrun
    
else
    echo "ERROR: File constant does not exit"
    exit -1
fi    

