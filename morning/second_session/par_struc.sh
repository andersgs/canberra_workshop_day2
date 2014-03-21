#!/bin/sh

START=$(date +%s)

maxRep=$1

task_id=$2

K=$[1+($task_id/$maxRep)]

rep=$[($K*$maxRep)-$task_id]

out=`echo output_K"$K"_r"$rep"`

chain=`echo chain_K"$K"_r"$rep"`

seed=`eval od -vAn -N4 -tu4 < /dev/urandom`

structure -K $K -o $out -D $seed > $chain

END=$(date +%s)

DIFF=$(echo "$END - $START" | bc)

echo "Finished job K=$K, replicate=$rep, and took $DIFF seconds. It had seed $seed"

