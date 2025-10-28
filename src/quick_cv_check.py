#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'

"""
This script is used to check if the cv devlops as the time goes by.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from ploting_set import acs_plot_style


def set_parser():
    parser = argparse.ArgumentParser(description="Check CV development over time")
    ### Parser stuff ###
    parser = argparse.ArgumentParser(description='calculate the free energy surfase (FES) along the chosen collective variables (1 or 2) using a reweighted kernel density estimate')
    # files
    parser.add_argument('-f',dest='filename',type=str,default='COLVAR',help='the COLVAR file name, with the collective variables and the bias')
    parser.add_argument('-o',dest='outfile',type=str,default='fes-rew.dat',help='name of the output file')
    # compulsory
    #parser.add_argument('-s',dest='sigma',type=str,required=True,help='the bandwidth for the kernel density estimation. Use e.g. the last value of sigma from an OPES_METAD simulation')
    # input columns
    parser.add_argument('--cv',dest='cv',type=str,default='2',help='the CVs to be used. Either by name or by column number, starting from 1')
    parser.add_argument('--skiprows',dest='skiprows',type=int,default=0,help='skip this number of initial rows')
    return parser.parse_args()

def plot_cv(data, output):
    """
    Plot the CV development over time.
    """
    plt.figure(figsize=(10, 6))
    x = data[:, 0]
    y_size = data.shape[1] - 1

    for i in range(y_size):
        plt.plot(x, data[:, i + 1], label=f"CV {i + 1}")

    plt.xlabel("Time")
    plt.ylabel("CV Value")
    plt.title("CV Development Over Time")
    plt.legend()
    plt.savefig(output)
    plt.close()

def main():
    args = set_parser()
    print("")
    error = "Error: %s"


    dim=len(args.cv.split(','))
    if dim==1:
        dim2 = False
    elif dim==2:
        dim2 = True
    else:
        sys.exit(error%('only 1D and 2D are supported'))

    with open(args.filename,'r') as f:
        fields=f.readline().split()
        if fields[1]!='FIELDS':
            sys.exit(error%('no FIELDS found in "%s"'%args.filename))
        try:
            col_x=int(args.cv.split(',')[0])-1
            name_cv_x=fields[col_x+2]
        except ValueError:
            sys.exit(error%('cv "%s" not found'%args.cv.split(',')[0]))
    if dim2:
        try:
            col_y=int(args.cv.split(',')[1])-1
            name_cv_y=fields[col_y+2]
        except ValueError:
            sys.exit(error%('cv "%s" not found'%args.cv.split(',')[1]))
        
    print(' using cv "%s" found at column %d'%(name_cv_x,col_x+1))
    if dim2:
        print(' using cv "%s" found at column %d'%(name_cv_y,col_y+1))

    skip_rows = args.skiprows + 1
    
    # read file
    all_cols=[0, col_x]
    if dim2:
        all_cols=[0, col_x,col_y]

    all_cols.sort() #pandas iloc reads them ordered
    data=pd.read_table(args.filename, dtype=float, sep='\s+', comment='#', header=None, usecols=all_cols, skiprows=skip_rows)
    if data.isnull().values.any():
        sys.exit(error%('your COLVAR file contains NaNs. Check if last line is truncated'))

    cv_x=np.array(data.iloc[:,all_cols.index(col_x)])
    if dim2:
        cv_y=np.array(data.iloc[:,all_cols.index(col_y)])

    acs_plot_style(use_style=False)
    if dim2:
        fig, ax = plt.subplots(1, 3, figsize=(12, 3))
    else:
        fig, ax = plt.subplots(1, 2, figsize=(8, 3))

    ax[0].plot(data.iloc[:, 0], cv_x, label=name_cv_x)
    if dim2:
        ax[0].plot(data.iloc[:, 0], cv_y, label=name_cv_y)
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("CV Values")
    ax[0].set_title("CV Development Over Time")
    ax[0].legend()

    ax[1].hist(cv_x, bins=101, density=True, alpha=0.5, label=name_cv_x)
    if dim2:
        ax[1].hist(cv_y, bins=101, density=True, alpha=0.5, label=name_cv_y)
    ax[1].set_xlabel("CV values")
    ax[1].set_ylabel("Density")
    ax[1].set_title("CV Distribution")
    ax[1].legend()

    if dim2:
        ax[2].scatter(cv_x, cv_y, s=1, alpha=0.5)
        ax[2].set_xlabel(name_cv_x)
        ax[2].set_ylabel(name_cv_y)
        ax[2].set_title("CV Scatter Plot")
    plt.savefig(args.outfile, bbox_inches='tight')

if __name__ == '__main__':
    main()

