#!/bin/bash

# Perform a full range of bench. First one much less intensive than second one
# which is closer to our real-application needs.

echo "Excitation frequency: 10, Sample rate: 200"
for i in 1 3 5 7 9 10 11 15 20 25 30
do
   echo $i
   ./RK4_LargeSystem.out -latticesize $i -ffrequency 10 -fsampling 200 -benchmarking | grep 'time'
done

echo "Excitation frequency: 72000, Sample rate: 4"
for i in 1 2 3 5 7
do
   echo $i
   ./RK4_LargeSystem.out -latticesize $i -ffrequency 72000 -fsampling 4 -benchmarking | grep 'time'
done