#!/bin/bash


mkdir out
NUM_MEMBERS=1
MIN_N=100
MAX_N=1000
MIN_SEED=1
MAX_SEED=10000

for i in {0..10..2}
  do
     echo "Welcome $i times"
 done

for numel in {0..100..10} 
do 
    echo $numel
done

for numel in {0..100..10} 
do 
    for seed in {1..10000..1} 
    do 
        for angle in {0..270..90} 
        do 
            for laneWidth in -1 2 4 
            do 
            xmlName=seed${seed}_angle${angle}_laneWidth${laneWidth}_numPoints${n}_numMembers${NUM_MEMBERS}
            echo ruby generate_random_weather.rb -seed ${seed} -lanewidth ${laneWidth} -angle ${angle} -operflex .82 2 4 -num_weather_points n -num_members ${NUM_MEMBERS} -oname out/${xmlName}
            echo RobustTree.exe -cinput inputs.txt
echo $numel
            done
        done
    done
done
