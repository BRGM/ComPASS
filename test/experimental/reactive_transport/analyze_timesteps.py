#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:21:22 2019

@author: kern
"""

import numpy as np
import matplotlib.pyplot as plt

def process(filename):
    dts = []
    for line in open(filename):
        if line[0:7] == 'Success':
            dts.append(line.split()[3])    
    dts = np.array(dts,dtype=np.double)
    ts = np.cumsum(dts)
    return ts, dts
        
if __name__=='__main__':
    process('momas_transp.out')
    print(ts, dts)
    plt.plot(ts, dts, '+-')
