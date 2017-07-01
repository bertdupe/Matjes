# this script makes the link between the movie.py script and the function that makes the spin
python movie.py
# ffmpeg is a very annoying way to encode (but it is the only one I know)
# the problem comes with -c:v libx264
# this part makes the encoding. To have a list of encoding ffmpeg -codecs
# if you have control over the computer, you can encode it as a codec libx264 or mp4 or whatever suits you
# the most commom encoding that I have found so far is png donc -c:v png
ffmpeg -r 10 -i spin_%d.png  -c:v libx264 -vf fps=12 -pix_fmt yuv420p anim.mp4
rm -rf *png
