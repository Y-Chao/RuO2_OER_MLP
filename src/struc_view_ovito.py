#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = 'Chao Yang'
__version__=	'1.0'

"""
Render a structure using OVITO and Tachyon.
"""

import os
import argparse
from typing import TYPE_CHECKING
from ase.io import read
from ovito.io import import_file
from ovito.modifiers import (CreateBondsModifier,
                             ReplicateModifier,
                             AffineTransformationModifier,
                             WrapPeriodicImagesModifier
                            )
from ovito.vis import Viewport, TachyonRenderer

if TYPE_CHECKING:
    from ovito.data import DataCollection

def modify_one_frame(frame: int, data: DataCollection):
    # Modify the data for the current frame
    
    # types list[int], where the length equals to the number of particles,
    #                        the value is the type of the particle.
    types = data.particles_.particle_types_

    for i in range(len(types)):
        if types[i] == 1:
            types.type_by_id_(types[i]).radius = 0.68 # 1 for Ru
        elif types[i] == 2:
            types.type_by_id_(types[i]).radius = 0.35 # 2 for O
        elif types[i] == 3:
            types.type_by_id_(types[i]).radius = 0.25 # 3 for H
            types.type_by_id_(types[i]).color = (0.9, 0.74, 0.74)

    # Display the cell boundary or not 
    data.cell_.vis.enabled = False

def create_bonds():
    # To do: modify the input arguments to a dictionary
    # Uniform, VdWRadius, Pairwise
    cb = CreateBondsModifier(mode =CreateBondsModifier.Mode.Pairwise)
    cb.vis.width = 0.2
    cb.set_pairwise_cutoff("H", "O", 1.1)
    cb.set_pairwise_cutoff("O", "Ru", 2.2)
    return cb

def set_viewport(viewer: str = "front"):
    # Set the viewport
    vp = Viewport()
    if viewer == "front":
        vp.type = Viewport.Type.Front
    elif viewer == "back":
        vp.type = Viewport.Type.Back
    elif viewer == "left":
        vp.type = Viewport.Type.Left
    elif viewer == "right":
        vp.type = Viewport.Type.Right
    elif viewer == "top":
        vp.type = Viewport.Type.Top
    elif viewer == "bottom":
        vp.type = Viewport.Type.Bottom
    elif viewer == "perspective":
        vp.type = Viewport.Type.Perspective
    elif viewer == "ortho":
        vp.type = Viewport.Type.Ortho
    else:
        raise ValueError("Invalid viewer type. Choose from: front, back, left, right, top, bottom, perspective, ortho.")
    vp.zoom_all()
    return vp

def parse_args():
    parser = argparse.ArgumentParser(description='Render a structure using OVITO and Tachyon.')
    parser.add_argument('-i', '--input', type=str, help='Input path to the structure file')
    parser.add_argument('-s', '--start', type=int, default=0, help='Start frame number')
    parser.add_argument('-e', '--end', type=int, default=0, help='End frame number')
    parser.add_argument('-f', '--frame', type=int, default=0, help='Frame number to render')
    parser.add_argument('--shift', type=float, default=0.0, help='Shift the structure')
    parser.add_argument('-o', '--output', type=str, default='output', help='Output image file name')
    return parser.parse_args()

def main():
    args = parse_args()
    filename = args.input
    start = args.start
    end = args.end
    shift_z = args.shift
    frame = args.frame
    output = args.output

    # Load the structure file
    pipeline = import_file(filename)
    # Modify the data for the current frame
    pipeline.modifiers.append(modify_one_frame)
    # Create bonds
    pipeline.modifiers.append(create_bonds())
    # Create repeated structures
    # pipeline.modifiers.append(ReplicateModifier(num_x=2, num_y=1, num_z=1))
    # Set the shift
    shift_z = AffineTransformationModifier(
        operate_on = {'particles'},
        transformation=[[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 1, shift_z]]
    )
    pipeline.modifiers.append(shift_z)
    # Wrap the structure
    pipeline.modifiers.append(WrapPeriodicImagesModifier())

    # Set the viewport
    pipeline.add_to_scene()
    vp = set_viewport(viewer="front")

    # Set the renderer
    renderer = TachyonRenderer(
        shadows = False,
        direct_light_intensity = 1.5,
        antialiasing = True,
        ambient_occlusion = True,
        ambient_occlusion_brightness = 0.9,
        ambient_occlusion_samples = 64,
        antialiasing_samples= 16
    )


    for i in range(start, end, frame):
        vp.render_image(size=(3840, 2160), filename=f"{output}-{i}.png", background=(1, 1, 1), renderer=renderer, frame=i)

if __name__ == '__main__':
    main()
