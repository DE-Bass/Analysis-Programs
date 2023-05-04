#!/bin/bash

echo "Starting SNID process"

for file in $(ls brad)
do
    if [[ $file == *.dat ]] # Get only .dat files
    then
        echo "Running SNID for:" $file

        # Move to info directory
        cd info
        
        # Run SNID for the file, no redshit adjustment, immediately quit after displaying
        export X="../brad/$file"
        printf 'n\nq' | /pkg/linux/snid/bin/snid zmin=0 zmax=0.15 plot=0 $X # Replace ../test_command.sh with /pkg/linux/snid/bin/snid

        # Exit into main directory
        cd ..
    fi
done

