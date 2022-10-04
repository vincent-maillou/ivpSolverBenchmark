#!/bin/bash

# Perform a full range of bench. First one much less intensive than second one
# which is closer to our real-application needs.

echo "Excitation frequency: 1, Sample rate: 100, Duration: 30"
for i in 1 3 5 7 9 10 11 15 20 25 30
do
   echo $i
   ./RK4_LargeSystem.out -latticesize $i -ffrequency 1 -pointsPerPeriode 100 -simDuration 30 -benchmarking | grep 'time'
done

echo "Excitation frequency: 72000, Sample rate: 4"
for i in 1 2 3 5 7
do
   echo $i
   ./RK4_LargeSystem.out -latticesize $i -ffrequency 72000 -pointsPerPeriode 4 -simDuration 1 -benchmarking | grep 'time'
done