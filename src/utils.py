#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'


import numpy as np

def exponential_moving_average(data, alpha=0.1):
    """Exponential moving average"""
    ema = np.zeros_like(data)
    ema[0] = data[0]
    for i in range(1, len(data)):
        ema[i] = alpha * data[i] + (1 - alpha) * ema[i - 1]
    return ema

def move_average(data, window_size):
    """
    Calculate the moving average of a given data.
    """
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')

def weighted_moving_average(data, weights):
    """
    Calculate the weighted moving average of a given data.
    """
    weights = np.array(weights)
    norm_weights = weights / np.sum(weights)
    return np.convolve(data, norm_weights, mode='same')

def gaussian_weights(window_size, sigma=None):
    if sigma is None:
        sigma = window_size / 6  # Default value, the windows range is +-3 sigma
    half = window_size // 2
    x = np.arange(-half, half + 1)
    weights = np.exp(-x**2 / (2 * sigma**2))
    return weights