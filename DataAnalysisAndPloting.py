#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 14:09:54 2022

@author: vincent-maillou
"""

import csv
import matplotlib.pyplot as plt
import subprocess
import math


""" Executing the simulation parsing the chosen parameter """

latticeSize = 2

sampleRate = 1000 # [Hz] reel
frequencyOfExcitationSource = 0.1 # [Hz] reel
numberOfPeriodeToSimulate = 1 # [N] reel

def NumberOfStepsToPerfom(sampleRate, numberOfPeriodeToSimulate, frequencyOfExcitationSource):
    return int(numberOfPeriodeToSimulate*(sampleRate/frequencyOfExcitationSource))

numberOfSteps = NumberOfStepsToPerfom(sampleRate, numberOfPeriodeToSimulate, frequencyOfExcitationSource) # [N] uint


bashCommand = "make"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = "./RK4_LargeSystem.out -out -latticesize " + str(latticeSize) + " -sampleRate " + str(sampleRate) + " -numberOfSteps "  + str(numberOfSteps) + ' -frequencyOfExcitationSource ' + str(frequencyOfExcitationSource)
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()



"""" Reading data from the output of the simulation and ploting them """

SIM_F = []
with open('SIM_F.csv', 'r') as SIM_F_file:
    SIM_F_csv = csv.reader(SIM_F_file)
    SIM_F_header = next(SIM_F_csv)
    
    for row in SIM_F_csv:
        SIM_F.append(float(row[0]))

SIM_X = []
with open('SIM_X.csv', 'r') as SIM_X_file:
    SIM_X_csv = csv.reader(SIM_X_file)
    
    SIM_X_header = next(SIM_X_csv)
    
    for i in range(len(SIM_X_header)-1):
        SIM_X.append([])
    
    for row in SIM_X_csv:
        for i in range(len(row)-1):
            SIM_X[i].append(float(row[i]))

SIM_F_t = []
SIM_F_t = [(i/len(SIM_F))*numberOfSteps/sampleRate for i in range(len(SIM_F))]

SIM_X_t = []
SIM_X_t = [(i/len(SIM_X[0]))*numberOfSteps/sampleRate for i in range(len(SIM_X[0]))]

# SIM_test = [0.36 for i in range(len(SIM_X[0]))]
# SIM_test_t = [(i/len(SIM_X[0]))*simDuration for i in range(len(SIM_X[0]))]

plt.figure()
for i in range(len(SIM_X)+1):
    if i == 0:
        plt.subplot(len(SIM_X)+1, 1, i+1)
        plt.ylabel('F(t) = cos(f(t))')
        plt.plot(SIM_F_t, SIM_F, 'r')
    else:
        plt.subplot(len(SIM_X)+1, 1, i+1)
        plt.ylabel('X_' + str(i) +'(t)')
        plt.plot(SIM_X_t, SIM_X[i-1])
        # plt.plot(SIM_test_t, SIM_test)
        
plt.xlabel('t [s]')
plt.show()

print("Max of X[0]: " + str(max(SIM_X[0])))
print("Min of X[0]: " + str(min(SIM_X[0])))
print("Max of F(t): " + str(max(SIM_F)))

    
    
    
    