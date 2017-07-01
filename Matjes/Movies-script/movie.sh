#!/bin/bash

python singlepanel.py

## libx264 encoder
ffmpeg -r 10 -i spin_%6d.png  -c:v libx264 -vf fps=25 -pix_fmt yuv420p anim-$data.mp4

## mpeg encoder
#ffmpeg -r 1 -i spin_%6d.png  -vcodec mpeg4 animation.avi

