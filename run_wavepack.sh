#!/bin/bash

# THIS SCRIPT WILL RUN A NEW WAVEPACK DYNAMICS SIMULATION. THE SIMULATION NEEDS AN INPUT FILE TO BE ABLE TO RUN. FIRST, THIS SCRIPT ASKS THE USER FOR THE PATH TO THE INPUT FILE CONTAINING ALL THE SIMULATION PARAMETERS. IF THE FILE DOES NOT EXISTS, THEN THE SCRIPT WILL GENERATE A NEW FILE USING THE INFORMATION GIVEN BY THE USER.


if [[ -z $1 ]]
then
   echo "You cannot run Wavepack_1D without an input file. Please give the input file as an argument."
   echo "If no input file exist, the argument should indicate the location where the file must be created"
   exit
else
   file=$1
   if [[ -f ${file} ]]
   then
      echo "Run Wavepack_1D with the existing input file ${file} ? (y or n)"
      read answer
      if [[ "${answer}" = 'y' || "${answer}" = 'n' ]]
      then
#         echo "Running Wavepack_1D or not... That is the question!"
      else
        while [[ "${answer}" != 'y' && "${answer}" != 'n' ]]
        do
           echo "Please answer by y or n. Do you want to run Wavepack_1D with the existing input file ${file} ?"
           read answer
        done
      fi
      if [[ "${answer}" = 'y' ]]
      then
         echo "Running the code!"
         #run the code... in principle
         exit
      fi
      if [[ "${answer}" = 'n' ]]
      then
         echo "Generation of a new input file..."
      fi
   else
      echo "Generation of a new input file..."
   fi
fi

############################
#THIS PART OF THE SCRIPT GENERATES THE NEW INPUT FILE
############################

echo "Accessing the second part of the script!"


exit
