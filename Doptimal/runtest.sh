#!/bin/bash
timelimit=3600
gnuparalleltest=1 # 1: use GNU parallel to speed up test; 0: not use
algorithms=("default" "icutss" "icutsl" "icutssa" "icutsla")
datapath="benchmark"
logpath="logs"
settingpath="settings"

# create and clear test files
#python3 ./testdir.py

runInstance() {
    timelimit=$1
    instance=$2
    algo=$3
    logpath=$4
    settingpath=$5
    datapath=$6

    #echo $timelimit $instance $benchmark $benchmarkpath $logpath $resultpath $algo
    

    #echo $logpath $instance $algo $instance $algo
    solver/build/dopt -c "set limits time $timelimit" -c  "set load $settingpath/$algo.set" -c  "read $datapath/$instance" -c "opt write statistics" -c "$logpath/${instance}_${algo}.log"  -c "quit"

}
export -f runInstance


instances=$(ls $datapath/$benchmark)
if [ $gnuparalleltest == 0 ]
then
    for instance in  $instances
    do
        for algo in ${algorithms[@]}
            do
                runInstance "$timelimit" "$instance" "$algo" "$logpath" "$settingpath" "$datapath"
        done
    done
else
    parallel --will-cite --jobs 85% runInstance  "$timelimit" ::: "$instances" :::  "${algorithms[@]}" ::: "$logpath" ::: "$settingpath" ::: "$datapath"
fi



