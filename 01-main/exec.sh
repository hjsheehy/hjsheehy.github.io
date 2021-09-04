#!/usr/bin/bash
confFolder=../02-conf/
dataFolder=../03-data/
plotFolder=../04-plot/
figFolder=../05-fig/
activeFolder=active/
inFolder=in/
outFolder=out/

num=0 # initialise var
pwr=0 
while getopts n:p: flag
do
    case "${flag}" in
        n) num=${OPTARG};;
        p) pwr=${OPTARG};;
    esac
done

if [ "$num" = "0" ]
then
    echo "Please supply simulation number (format: -n ##-##)"
    exit 2
fi
model="${num:0:2}"
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
echo    # (optional) move to a new line
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
echo -e "${confOutFolder}\nConfiguration files found"

    read -p 'Delete and create new .conf files? [Y/N]' -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        python3 $confFolder/conf.py -s
    fi
    
fi    

if [ `ls $confOutFolder | wc -l` -eq 0 ] 
then
    echo -e "${confOutFolder}\nConfiguration out folder empty"
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
    echo -e "${dataFolder}\nFound data (.npz files)."
    read -p 'Delete? [Y/N]' -n 1 -r
    echo    # (optional) move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        rm -r ${dataFolder}
        mkdir $dataFolder
    fi
fi    

plt=false
if [ ${pwr} = '0' ]
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
fi

if ${plt}; then
    echo `xdg-open $figFolder` &
fi
function func {
    while [ `ls $inFolder | wc -l` -gt 0 ] ; do
        grab_name=`ls ${inFolder}/*.conf | head -n 1`
        grab=`basename ${grab_name}`
        mv ${inFolder}/${grab} ${activeFolder}/${grab}
        echo "Running main.py ${grab}"
        python3 $model*.py ${activeFolder}/${grab}
        wait
        echo "Done main.py ${grab}"

        mv ${activeFolder}/${grab} ${confOutFolder}

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
