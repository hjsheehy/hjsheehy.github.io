#! /bin/bash

Main(){
  i=1
  while [ $i -le $n_nodes ] 
  do
    qsub self-consistent.pbs
    sleep 10
   ((i++))
  done
  echo ALL done
}

n_nodes=13
Main