#!/usr/bin/bash

# PythU Unix shell scripting module
# September 12th, 2021
__version__='1.0.0'

# A module for scripting large PythU simulations
#
# Copyright (C) 2021, Henry Joseph Sheehy
# 
# PythU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PythU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# A copy of the GNU General Public License should be available
# alongside this source in a file named TBC.  If not,
# see <http://www.gnu.org/licenses/>.
#
# PythU is availabe at http://www. TBC

####################################################################
# License                                                          #
####################################################################
License()
{
    # Display License
    echo "PythU Unix shell scripting module"
    echo "September 12th, 2021"
    echo "__version__='1.0.0'"
    echo 
    echo "A module for scripting large PythU simulations"
    echo 
    echo "Copyright (C) 2021, Henry Joseph Sheehy"
    echo 
    echo "PythU is free software: you can redistribute it and/or modify"
    echo "it under the terms of the GNU General Public License as published by"
    echo "the Free Software Foundation, either version 3 of the License, or"
    echo "(at your option) any later version."
    echo 
    echo "PythU is distributed in the hope that it will be useful,"
    echo "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    echo "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    echo "GNU General Public License for more details."
    echo 
    echo "A copy of the GNU General Public License should be available"
    echo "alongside this source in a file named TBC.  If not,"
    echo "see <http://www.gnu.org/licenses/>."
    echo 
    echo "PythU is availabe at http://www. TBC"
    echo
}
####################################################################
# Help                                                             #
####################################################################
Help()
{
   # Display Help
   echo "Add description of the script functions here."
   echo
   echo "Syntax: scriptTemplate [-g|h|n|p]"
   echo "options:"
   echo "g     Print the GPL license notification."
   echo "h     Print this Help."
   echo "n     Select simulation number (format: ##-##)"
   echo "p     Optional power option after programme ended."
   echo "      options: [shutdown|sleep]"
   echo
}

num="" # initialise var
pwr="" 

while getopts :hgn:p: flag; do
   case ${flag} in
      h) # display Help
         Help
         exit;;
      g) # display License
         License
         exit;;
      n) num=${OPTARG};;
      p) pwr=${OPTARG};;
     \?) echo "Error: Invalid option"
         exit;;
    esac
done
####################################################################
# Main                                                             #
####################################################################
confFolder=../02-conf/
dataFolder=../03-data/
plotFolder=../04-plot/
figFolder=../05-fig/
activeFolder=active/
inFolder=in/
outFolder=out/

# while getopts n:p flag
# do
#     case "${flag}" in
#         n) num=${OPTARG};;
#         p) pwr=${OPTARG};;
#     esac
# done

if [ "$num" = "" ]
then
    echo "Please supply simulation number (format: -n ##-##) or call -h for help."
    exit 2
fi
model="${num:0:2}"
model=($model*.py)
#########################################################
#################### power options ######################
#########################################################
sdn=false
slp=false

if [ $pwr = "sleep" ]; then
    echo "Sleeping after completion."
    slp=true
fi

if [ $pwr = "shutdown" ]; then
    echo "Shutting down after completion."
    sdn=true
fi
##########################################################
###################### architecture ######################
##########################################################
filename=$(basename $confFolder$num*)

confFolder=$confFolder/$filename
confOutFolder=$confFolder/$outFolder
dataFolder=$dataFolder/$filename
plotFolder=$plotFolder/$filename
figFolder=$figFolder/$filename.pdf
activeFolder=$activeFolder/$filename
inFolder=$inFolder/$filename
# outFolder=$outFolder/$filename

echo 'Simulation: '$filename
read -p 'Continue? [Y/N]' -n 1 -r
if [[ $REPLY =~ ^[Nn]$ ]]
then
    exit 1
fi

mkdir -p $activeFolder $inFolder $confOutFolder

if [ `ls $activeFolder | wc -l` -gt 0 ] 
then
    echo `rm -r ${activeFolder}`
fi
if [ `ls $inFolder | wc -l` -gt 0 ] 
then
    echo `rm -r ${inFolder}`
fi
# if [ `ls $outFolder | wc -l` -gt 0 ] 
# then
#     echo `rm -r ${outFolder}`
# fi

mkdir -p $confFolder $dataFolder $activeFolder $inFolder 

if [ `ls $confOutFolder | wc -l` -gt 0 ] 
then
echo -e "Configuration files found in\n${confOutFolder}"

    read -p 'Delete and create new .conf files? [Y/N]' -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        python3 $confFolder/conf.py -s
    fi
    
fi    

if [ `ls $confOutFolder | wc -l` -eq 0 ] 
then
    echo -e "Configuration out folder empty, located at\n${confOutFolder}"
    read -p 'Create .conf files? [Y/N]' -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Nn]$ ]]
    then
        echo 'No configuration files. Exiting...'
        exit 1
    else
    python3 $confFolder/conf.py -s
    fi
fi    

if [ `ls $dataFolder | wc -l` -gt 0 ] 
then
    echo -e "Found data (.npz files) at\n${dataFolder}"
    read -p 'Delete? [Y/N]' -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        rm -r ${dataFolder}
        mkdir $dataFolder
    else
        for f in ${dataFolder}/*.npz ; do
            f="${confOutFolder}/`basename ${f%.*}.conf`"
            mv $f $confFolder
        done
    wait
    fi
fi    

plt=false
if [ ${pwr} = '' ]
then
    read -p 'Plot output? [Y/N]' -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        plt=true
    fi
fi

########################################################
##################### scripting ########################
########################################################
echo 'Execution beginning...'

if [ `ls $confOutFolder | wc -l` -gt 0 ] 
then
  mv $confOutFolder/* $inFolder
  mv $confFolder/*.conf $confOutFolder
fi

if ${plt}; then
    echo `xdg-open $figFolder` &
fi
function func {
    while [ `ls $inFolder | wc -l` -gt 0 ] ; do
        grab_name=`ls ${inFolder}/*.conf | head -n 1`
        grab=`basename ${grab_name}`
        mv ${inFolder}/${grab} ${activeFolder}/${grab}
        wait
        echo "Running ${model} ${grab}"
        python3 ${model} ${activeFolder}/${grab}
        wait
        echo "Done ${model} ${grab}"

        mv ${activeFolder}/${grab} ${confOutFolder}
        wait

        if ${plt}
        then
            python3 ${plotFolder}.py -s > /dev/null 2>&1 &
        fi
    done
}
  
function main {
    for ((i=0 ; i<$n_parallel ; i++)); do
        func &
    done
    wait
    rm -r $inFolder $activeFolder
    echo "Execution finished"

    if ${slp}
    then
        echo `systemctl suspend`
    fi

    if ${sdn}
    then
        echo `systemctl poweroff`
    fi
}

#######################################################
####################### Main ##########################
#######################################################
n_parallel=1
main
