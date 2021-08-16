#!/bin/bash

path=$(pwd)
echo "Hello there! You are using the MATJES Code (let's be honest: this is awesome)"
echo "This script should will try to create a standard top image out of your calculated magnetic spin structure"
echo "If you do not wish run it automatically, adapt the povrayscript.pov file and run the command"
echo "-->"
echo "povray -w1000 -h1000 povrayscript.pov"
echo "Here, -w1000 -h1000 defines the resolution for width and height"
echo ""

loc_x=$(cat positions.dat | gawk '{print $1}' | LC_ALL=C sort -k 1 -g | tail -n 1 | gawk '{print -$1/2}')
loc_y=$(cat positions.dat | gawk '{print $2}' | LC_ALL=C sort -k 1 -g | tail -n 1 | gawk '{print $1/2}')
totalheight=$(cat positions.dat | gawk '{print $3}' | LC_ALL=C sort -k 1 -g | tail -n 1)
loc_z=$(echo "$loc_x  $loc_y  $totalheight" | gawk '{print 1.5*sqrt($1*$1+$2*$2+$3*$3)}')

if [ ! -d "Visualization" ]
then
mkdir -p $path/Visualization
fi
###########################
# prepare file for povray #
###########################

for state in start end
do
paste positions.dat magnetic_$state.dat | gawk '{for (i=1;i<=NF;i++) printf "%.10f%s", $i, (i<NF?OFS:ORS)}' | gawk '{print "Spin(", "\t", $1",\t",  $2",\t", $3",\t", $4",\t", $5",\t", $6")"}' > Visualization/position_magnetization_$state.dat
done


cp povrayscript.pov $path/Visualization
cd $path/Visualization

sed -i "/camera {/{n;s/.*/location < $loc_x, 0, $loc_z >\nlook_at  < $loc_x, 0, 0 > /}" povrayscript.pov

for state in start end
do
sed -i "$ s/#include .*/#include \"position_magnetization_$state.dat\"/g" povrayscript.pov
##############
# run povray #
##############
povray -w1000 -h1000 povrayscript.pov Output_File_Name=spinstructure_$state.png

done
echo "Thank you very much!!!"
echo "Your pictures are now in the directory 'Visualization'"
