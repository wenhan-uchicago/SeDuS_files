#!/bin/bash

# B is the number of blocks
B=3;
# name is the name identifier used in SeDuS
name="ha";

# the (SUPERTIME / multi_thread_limit) should be an integer
multi_thread_limit=10;

SUPERTIME=10;

count_time=1;

while [ $SUPERTIME -gt 0 ]
do
    i=1;

    # running a number of sub-processes
    while [ $i -le $multi_thread_limit ]
    do
	echo "Real supertime $count_time $i";
	echo "Real supertime $count_time $i" >> temp_output_${i};
	./SeDuS_original ${name}_${i} >> temp_output_${i} &

	# necessary to sleep a random amount of time; otherwise results may be the same
	sleep 1;
	# vital, otherwise server will die
	let i=$i+1;

    done

#    sleep 5;

    i=1;

    while [ $i -le $multi_thread_limit ]
    do
	
	# test if all sub-processes have finished
	judge_finished=`tail -n 1 temp_output_${i}`;

	# "" vital!
	while [ "${judge_finished}" != "Finished Run" ]
	do
	    echo "$i Not finished yet; $judge_finished";
	    sleep 5;
	    judge_finished=`tail -n 1 temp_output_${i}`;
	done
	echo "$i finished already; $judge_finished";

	# if finished, merge all information	
	j=0;
	while [ $j -lt $B ]
	do
	    cat samplepi${j}[50]_${name}_${i}.dat >> total_samplepi${j}[50]_${name}.dat;
	    cat sampleS${j}[50]_${name}_${i}.dat >> total_sampleS${j}[50]_${name}.dat;
	    cat mutations_new_${j}_${name}_${i}.dat >> total_mutations_new_${j}_${name}.dat;
	    cat SFS_${j}_${name}_${i}.dat >> total_SFS_${j}_${name}.dat;

	    # remove files here
	    #rm samplepi${j}[50]_${name}_${i}.dat;
	    #rm sampleS${j}[50]_${name}_${i}.dat;
	    #rm mutations_new_${j}_${name}_${i}.dat;
	    #rm SFS_${j}_${name}_${i}.dat;

	    let j=$j+1;
	done

	# merge all output files
	cat temp_output_${i} >> multithread_output;
	
	# remove those files


	let i=$i+1;
	let count_time=$count_time+1;
    done

    let SUPERTIME=$SUPERTIME-$multi_thread_limit;
done
