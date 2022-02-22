#!/usr/bin/bash

# PythU Unix shell scripting module
# January, 2022
__version__='1.0.0'

# A module for scripting large PythU simulations
#
# Copyright (C) 2022, Henry Joseph Sheehy
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
    echo "January, 2022"
    echo "__version__='1.0.0'"
    echo 
    echo "A module for scripting large PythU simulations"
    echo 
    echo "Copyright (C) 2022, Henry Joseph Sheehy"
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
   echo "Syntax: scriptTemplate [-g|h|s|n|c|p]"
   echo "options:"
   echo "g     Print the GPL license notification."
   echo "h     Print this Help."
   echo "s     Select simulation number (format: ##-##)"
   echo "n     Select number of parallel simulations to run. (If c=true, then each simulation is run on a seperate node on the cluster.)"
   echo "c     supercomputing cluster (PBS)?" 
   echo "      options: [true|false]"
   echo "p     Optional power option after programme ended."
   echo "      options: [shutdown|sleep]"
   echo
}

clr=false
num="" # initialise var
pwr="init" 
nodes=1

while getopts :hgs:n:c:p: flag; do
   case ${flag} in
      h) # display Help
         Help
         exit;;
      g) # display License
         License
         exit;;
      s) num=${OPTARG};;
      n) nodes=${OPTARG};;
      c) clr=${OPTARG};;
      p) pwr=${OPTARG};;
     \?) echo "Error: Invalid option"
         exit;;
    esac
done

if $clr; then
    module load python3
fi

function Main {

	if [ "$num" = "" ]
	then
	    echo "Please supply simulation number (format: -n ##-##) or call -h for help."
	    exit 2
	fi
	modelConf="${num:0:5}"
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
	plotFolder=$plotFolder$modelConf
	figFolder=$figFolder$modelConf
	activeFolder=$activeFolder/$filename
	inFolder=$inFolder/$filename
	# outFolder=$outFolder/$filename

	echo 'Simulation: '$filename
	read -p 'Continue? [Y/N]' -n 1 -r
	if [[ $REPLY =~ ^[Nn]$ ]]
	then
	    exit 1
	fi

	mkdir -p $activeFolder $inFolder $confFolder $dataFolder 
	mkdir -p $confOutFolder 


	if [ `ls $inFolder | wc -l` -gt 0 ] 
	then
	    if [ `ls $activeFolder | wc -l` -gt 0 ]; then
		echo -e "Detected simulation was stopped (data in \n${activeFolder})"
		read -p 'Continue simulation? [Y/N]' -n 1 -r
		echo
		if [[ $REPLY =~ ^[Yy]$ ]]
		then
		    mv $activeFolder/* $inFolder
		else
		    echo `rm -r ${activeFolder}/*`
		    echo `rm -r ${inFolder}/*`
		fi
	    fi
	fi

	if [ `ls $confOutFolder | wc -l` -gt 0 ] 
	then
	echo -e "Configuration files found in\n${confOutFolder}"

	    read -p 'Delete and create new .conf files? [Y/N]' -n 1 -r
	    echo    # (optional) move to a new line
	    if [[ $REPLY =~ ^[Yy]$ ]]
	    then
		${python} $confFolder/conf.py -s
	    else
		if [ `ls $activeFolder | wc -l` -gt 0 ] 
		then
		    echo `rm -r ${activeFolder}`
		fi
		if [ `ls $inFolder | wc -l` -gt 0 ] 
		then
		    echo `rm -r ${inFolder}`
		fi
		mkdir -p $confFolder $dataFolder $activeFolder $inFolder 
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
	        ${python} $confFolder/conf.py -s
	    fi
	fi    

	if [ `ls $dataFolder | wc -l` -gt 0 ] 
	then
	    echo -e "Found data (.npz files) at\n${dataFolder}"
	    read -p 'Delete? [Y/N]' -n 1 -r
	    echo    # (optional) move to a new line
	    if [[ $REPLY =~ ^[Yy]$ ]]
	    then
		echo `rm -r ${dataFolder}/*`
	    fi
	fi

	for fullfile in ${dataFolder}/*.npz
	do
	    filename=$(basename -- "$fullfile")
	    filename="${filename%_*}"
	    file=${confOutFolder}/$filename*
	    if [ -f $file ] ; then
		rm $file
	    fi
	done
	mv $confOutFolder/* $inFolder

	plt=false
	if [ ${pwr} = 'init' ]
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
	open_figure=${plt}

	function func {
	    while [ `ls $inFolder | wc -l` -gt 0 ] ; do
		grab_name=`ls ${inFolder}/*.conf | head -n 1`
		grab=`basename ${grab_name}`
		mv ${inFolder}/${grab} ${activeFolder}/${grab}
		wait
		echo "Running ${model} ${grab}"

		if $clr; then
		    bash qsub_python.sh ${python} ${model} ${activeFolder}/${grab}
		    while [ -f "${activeFolder}/${grab}" ] ; do
			sleep .$[ ( $RANDOM % 4 ) + 1 ]s
		    done
		    
		else
		    ${python} ${model} ${activeFolder}/${grab}
		    wait
		fi

		echo "Done ${model} ${grab}"

		wait

		if ${plt}
		then
		    for plot_script in ${plotFolder}*.py; do
			${python} $plot_script -s > /dev/null 2>&1 &
			if ${open_figure}; then
			    open_figure=false
			    (echo `xdg-open $figFolder*.pdf` &) &
			    echo $open_figure &
			fi
		    done
		fi
	    done
	}
	  
	    for ((i=0 ; i<${nodes} ; i++)); do
		func &
		sleep 2s
	    done
	    wait
	    rm -r $inFolder $activeFolder
	    echo "Execution finished"
	    notify-send "Execution finished"

	    if ! ${plt}
	    then
		for plot_script in ${plotFolder}*.py; do
		    ${python} $plot_script -s > /dev/null 2>&1 &
		done
	    fi

	    if ${slp}
	    then
		echo "Sleeping in 30 seconds..."
		echo "[ctrl+c] to cancel."
		notify-send "Sleeping in 30 seconds..."
		sleep 30s
		echo `systemctl suspend`
	    fi

	    if ${sdn}
	    then
		echo "Shuting down in 30 seconds..."
		echo "[ctrl+c] to cancel."
		notify-send "Shuting down in 30 seconds..."
		sleep 30s
		echo `systemctl poweroff`
	    fi
}
#######################################################
####################### Main ##########################
#######################################################
confFolder=../02-conf/
dataFolder=../03-data/
plotFolder=../04-plot/
figFolder=../05-fig/
activeFolder=active/
inFolder=in/
outFolder=out/
python=python3

Main
