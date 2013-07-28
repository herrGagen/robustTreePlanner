#!/bin/bash

mkdir out${1}-${2}
NUM_MEMBERS=3
MIN_N=100
MAX_N=1000
MIN_SEED=$1
MAX_SEED=$2

for i in {0..10..2}
  do
     echo "Welcome $i times"
 done

for numel in {0..100..10} 
do 
    echo $numel
done

for ((numEl = ${MIN_N}; numEl <= ${MAX_N}; numEl+= 100))
do 
    for ((seed  = ${MIN_SEED}; seed<=${MAX_SEED}; seed++))
    do 
        for angle in {0..270..90} 
        do 
            for laneWidth in -1 2 4 
            do 
            xmlName=seed${seed}_angle${angle}_laneWidth${laneWidth}_numPoints${numEl}_numMembers${NUM_MEMBERS}
            ruby generate_random_weather.rb -seed ${seed} -lanewidth ${laneWidth} -angle ${angle} -operflex .82 2 4 -num_weather_points ${numEl} -num_members ${NUM_MEMBERS} -oname out${1}-${2}/${xmlName}
            ./rtp -cinput inputs.txt
#echo $seed
            done
        done
    done
done
