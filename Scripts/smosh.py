#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:07:26 2020

@author: jacob
"""
import os

directory = '/home/jacob/examensarbete/software/output/5000/QC_Spiral(false)/dist(0.75)/5'
saveDirectory = '/home/jacob/examensarbete/Scripts/'
data = ''
for filename in os.listdir(directory):
    with open(os.path.join(directory, filename)) as file:
        for line in file.readlines():
            if line != 'gene_id\n': 
                data += line
                
with open(os.path.join(saveDirectory,"output.txt"), 'w') as file:
    file.write(data)