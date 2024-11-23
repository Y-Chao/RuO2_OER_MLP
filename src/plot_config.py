#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Chao Yang'
__version__=	'1.0'

import os
from scipy.constants import golden_ratio 
from matplotlib import pyplot as plt
from matplotlib import rcParams

def _plotting_config(width=3.37, height=3.37/golden_ratio*1.5/2,\
                     lrbt=[0.135,0.80,0.25,0.95], ffont='Arial',fsize=9.0):
    """
    Function to set parameters for formatting of standard plots

    Parameters
    ----------
    width: float
        Width of the plot in inches
    height: float
        Height of the plot in inches
    lrbt: list/tuple
        List of floats to set the left, right, bottom, and top margins
    ffont: string
        Font-type for the plot
    fsize: float
        Font-size for the plot
    """
    # set plot geometry
    rcParams['figure.figsize'] = width, height
    rcParams['figure.subplot.left'] = lrbt[0]
    rcParams['figure.subplot.right'] = lrbt[1]
    rcParams['figure.subplot.bottom'] = lrbt[2]
    rcParams['figure.subplot.top'] = lrbt[3]

    # set font properties
    rcParams['font.sans-serif'] = ffont
    rcParams['font.size'] = fsize

    # set tick properties   
    rcParams['axes.linewidth'] = 2
    rcParams['axes.labelsize'] = fsize

    # set tick properties
    rcParams['xtick.top'] = False # set to true to draw ticks on the top side
    #rcParams['xtick.direction'] = 'in'
    rcParams['ytick.right'] = False # set to true to draw ticks on the right side
    #rcParams['ytick.direction'] = 'in'
    rcParams['xtick.major.size'] = 8
    rcParams['xtick.major.width'] = 2
    rcParams['xtick.minor.size'] = 5
    rcParams['xtick.minor.width'] = 2
    rcParams['ytick.major.size'] = 8
    rcParams['ytick.major.width'] = 2
    rcParams['ytick.minor.size'] = 5
    rcParams['ytick.minor.width'] = 2
    rcParams['xtick.labelsize'] = fsize
    rcParams['ytick.labelsize'] = fsize

    # set legend properties
    rcParams['legend.fancybox'] = False
    rcParams['legend.edgecolor'] = 'k'

def writefig(filename, folder='output', write_pdf=False, write_eps=False, write_png=True):
    """
    Wrapper for creating figures, it should be used after plotting

    Parameters
    ----------
    filename: string
        Name of the file to save
    folder: string
        Name of the folder to save the file
    write_pdf: bool
        Flag to write a pdf file
    write_eps: bool
        Flag to write an eps file
    write_png: bool
        Flag to write a png file
    """
    # folder for output
    if not os.path.exists(folder):
        os.makedirs(folder)

    fileloc = os.path.join(folder, filename)
    if write_pdf:
        plt.savefig(fileloc+'.pdf')
    if write_eps:
        plt.savefig(fileloc+'.eps')
    if write_png:
        plt.savefig(fileloc+'.png', dpi=300)

    plt.close(plt.gcf())

clrs = {'orange':[0.8901960784313725, 0.4470588235294118, 0.13333333333333333],
        'darkgray':[0.34509803921568627, 0.34509803921568627, 0.35294117647058826],
        'lightblue':[0.596078431372549, 0.7764705882352941, 0.9176470588235294],
        'lightblue2':[0.39215686274509803, 0.6274509803921569, 0.7843137254901961],
        'lightblue3': [0.0, 0.6, 1.0],
        'deepblue':[0.0, 0.396078431372549, 0.7411764705882353],
        'azurblue':[0., 0.833, 1.],
        'darkblue':[0., 0., 0.5], 
        'byellow':[1., 0.90, 0.], 
        'darkyellow':[1.0, 0.706, 0.0],
        'orange2':[0.8980392156862745, 0.20392156862745098, 0.09411764705882353],
        'darkred':[0.5, 0., 0.],
        'green':[0.6352941176470588, 0.6784313725490196, 0.0],}
