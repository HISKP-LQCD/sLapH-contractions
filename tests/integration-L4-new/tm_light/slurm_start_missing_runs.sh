#!/bin/bash

if [ ! -e missing_configs.txt ]; then
  echo "missing_configs.txt (generated by check_perams.py) must exist in the pwd"
  kill -INT $$
else
  while IFS='' read -r line || [[ -n "$line" ]]; do
    # skip commented lines
    case $line in
      ''|\#*) continue;;
    esac

    cnfg=$( echo $line | awk -F ',' '{print $1}' )
    rnd=$( echo $line | awk -F ',' '{print $2}' )
    pushd cnfg$( printf %04d ${cnfg} )/rnd_vec_$(printf %02d $rnd)
      echo submitting cnfg=$cnfg rnd_vec=$rnd
      sbatch *job*.cmd
    popd
  done < missing_configs.txt
fi

