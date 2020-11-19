#!/bin/bash
W=$(grep resolution ../dyna.in | awk {'print$2'})
H=$(grep resolution ../dyna.in | awk {'print$3'})
rm -f *.mp4 toto.pov

for data in Ffield Spinse
do

nfile=$(ls $data* | grep -c $data)

if [ $nfile != 0 ]
then

end=$(($nfile+1))

name=$(grep include spin.pov | grep .dat | awk '{print$2}')

sed 's/'$name'/"'$data'_1.dat"/g' spin.pov > toto.pov

povray -W$W -H$H -V -UD -D +ospin_000001.png toto.pov 2> /dev/null
rm toto.pov

i=2
while ((i<end))
do
im=$((i-1))
sed 's/'$name'/"'$data'_'$i'.dat"/g' spin.pov > toto.pov
printf -v num "%06d" $i
povray -W$W -H$H -V -UD -D +ospin_"$num".png toto.pov 2> /dev/null
rm toto.pov
let i++
done
end=$(($end-1))
## libx264 encoder
ffmpeg -r 10 -i spin_%6d.png  -c:v libx264 -vf fps=25 -pix_fmt yuv420p anim-$data.mp4
rm spin_*

## mpeg encoder
#ffmpeg -r 1 -i spin_%6d.png  -vcodec mpeg4 animation.avi

fi
done
