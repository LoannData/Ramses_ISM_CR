#!/bin/bash
#MSUB -r tstranex6           # Nom du job                
#MSUB -n 16                  # Reservation de 4 processus au total
#MSUB -N 8                  # Les processus sont répartis sur 2 noeuds
#MSUB -T 36000                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_test5ramses          # Sortie standard
#MSUB -e stderr_test5ramses          # Sortie d'erreur       
#MSUB -p gen2191      # Allocation
#MSUB -q gpu                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x

export DATE=`date +%F_%Hh%M`

mpirun $BRIDGE_MSUB_PWD/ramses3d $BRIDGE_MSUB_PWD/test5_ramses.nml > log_$DATE.log