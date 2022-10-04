#!/bin/bash

# Perform a full range of bench. First one much less intensive than second one
# which is closer to our real-application needs.

echo "frequencyOfExcitationSource: 1, sampleRate: 1000, numberOfSteps: 10000"
for i in 1 3 5 7 9 10 11 15 20 25 30
do
   echo $i
   ./RK4_LargeSystem.out -latticesize $i -sampleRate 1000 -numberOfSteps 10000 -frequencyOfExcitationSource 1 -benchmarking | grep 'time'
done

echo "frequencyOfExcitationSource: 72000, sampleRate: 720000, numberOfSteps: 720000"
for i in 1 2 3 5 7
do
   echo $i
   ./RK4_LargeSystem.out -latticesize $i -sampleRate 720000 -numberOfSteps 720000 -frequencyOfExcitationSource 72000 -benchmarking | grep 'time'
done