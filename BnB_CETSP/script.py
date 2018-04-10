#!/usr/bin/env python
import os, sys, signal, fnmatch
from os import listdir
from os.path import isfile, join
from subprocess import call
from numpy import *
from random import shuffle

#print arguments
if(len(sys.argv) != 2):
    print 'Calling syntax: python script.py [PATH]';
    print 'Number of arguments:', len(sys.argv), 'arguments.';
    print 'Argument List:', str(sys.argv);
    exit(1);

#get path
myPath = sys.argv[1];

#get required matching files
allFiles = []
for file in listdir(myPath):
    if fnmatch.fnmatch(file, r'tp-[2][5]*') or fnmatch.fnmatch(file, r'tp-[3][0]*'):
        allFiles.append(file);

#sort list
allFiles.sort();
print allFiles

#define radii
radii = ['1.0', '0.5', '0.25'];

#main loop
for radius in radii:
    for instance in allFiles:
        print "Running instance: " + instance;
        command = "./exeCVXHULL Behdani/TPs/" + instance + " 2D " + radius + " 3600 V1 BeFS 1 1";
        print command;
        var = os.system(command);
        if (var == 2):
            print "\nKeyboard interruption!\n";
            sys.exit();
    os.system("mkdir Resultados/Wang/" + radius);
    os.system("mv Resultados/MATLAB/* Resultados/Wang/" + radius);
    os.system("mv Resultados/resultados_BnB_CETSP_OPTF.txt Resultados/Wang/" + radius);




