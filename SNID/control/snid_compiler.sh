#!/bin/bash

echo "Starting SNID process"

export Y=0

for file in $(ls brad)
do
    if [[ $file == *.dat ]] # Get only .dat files
    then
        # Construct the path of the info file
        info_file=${file::${#file}-4}
        info_file+="_snid.output"

        # Move to info directory
        cd info

        # Check that SNID does not already have an existing file
        if [ ! -e $info_file ]; then

            echo "Running SNID for:" $file
        
            # Run SNID for the file, no redshit adjustment, immediately quit after displaying
            export Y=$(($Y+1))
            echo $Y
            export X="../brad/$file"
            printf 'n\nq' | /pkg/linux/snid/bin/snid zmin=0 zmax=0.15 plot=0 $X # /pkg/linux/snid/bin/snid zmin=0 zmax=0.15 plot=0 $X # Replace ../test_command.sh with /pkg/linux/snid/bin/snid

            # Check that SNID saved the file
            if [ ! -e $info_file ]; then
                # If failed, rerun SNID with larger redshift
                printf 'n\nq' | /pkg/linux/snid/bin/snid zmin=0 zmax=1 plot=0 $X
            fi

            # Check again that SNID saved the file
            if [ ! -e $info_file ]; then
                # If failed again, lower rlapmin to zero
                printf 'n\nq' | /pkg/linux/snid/bin/snid zmin=0 zmax=1 rlapmin=0 plot=0 $X
            fi

        fi

        # Exit into main directory
        cd ..

    fi
done

# Check that all files were loaded
total1=$(ls brad -1 | grep "\.dat" | wc -l)
total2=$(ls info -1 | grep "\.output" | wc -l)
echo "Checksum: $Y files loaded, $total2/$total1 files saved."
