#!/bin/bash

path=$(pwd)
echo "Hello there! You are using the MATJES Code (let's be honest: this is awesome)"
echo "This script should will try to create a standard top image out of your calculated magnetic spin structure"
echo "If you do not wish run it automatically, adapt the povrayscript.pov file and run the command"
echo "-->"
echo "povray -w1000 -h1000 povrayscript.pov"
echo "Here, -w1000 -h1000 defines the resolution for width and height"
echo ""

# usage examples of -s|--states option:
# ./visualizationScript.sh -s "(start 25 50 75 end)"
# ./visualizationScript.sh -s "(start $(seq 25 25 75) end)"
# ./visualizationScript.sh -s "(start  end)" -os

xcut=false
cutpos_y=0
s_option=false
os=false
resolution=1000 # default resolution
while [ -n "$1" ]
do
    case "$1" in
    -x | --xcut) # cut along the x axis at y=cutpos_y=cste
        xcut=true
        if [[ -n "$2" ]]
        then cutpos_y=$2
        else
            echo "error not enough paramters given for -x/--xcut option"
        fi
        ;;
    -s | --states) # gives a list of states
        # if [[ -n "$2" ]]
        #     then echo "error not enough paramters given for -s/--states option"
        #     exit
        # fi
        s_option=true
        echo $2
        states=$2
        eval states=$2
        ;;
    -os) # output files have the input sates as suffixes
        os=true
        ;;
    -r | --resolution)
        resolution=$2
        echo "changed default resolution (1000x1000) to ${2}x${2}"
    esac
    shift
done

posx_min=$(cat positions.dat | gawk '{print $1}' | LC_ALL=C sort -k 1 -g | head -n 1)
posx_max=$(cat positions.dat | gawk '{print $1}' | LC_ALL=C sort -k 1 -g | tail -n 1)

posy_min=$(cat positions.dat | gawk '{print $2}' | LC_ALL=C sort -k 1 -g | head -n 1)
posy_max=$(cat positions.dat | gawk '{print $2}' | LC_ALL=C sort -k 1 -g | tail -n 1)

posz_min=$(cat positions.dat | gawk '{print $3}' | LC_ALL=C sort -k 1 -g | head -n 1)
posz_max=$(cat positions.dat | gawk '{print $3}' | LC_ALL=C sort -k 1 -g | tail -n 1)

if [ $xcut == true ]
then
loc_x=$(echo "$posx_min  $posx_max" | gawk '{print -($2-$1)/2}')
loc_y=$(echo "$loc_x  $cutpos_y" | gawk '{print -1.5*sqrt(2*$1*$1)+$2}')
loc_z=0
else
loc_x=$(echo "$posx_min $posx_max" | gawk '{print -($2-$1)/2}')
loc_y=$(echo "$posy_min $posy_max" | gawk '{print ($2-$1)/2}')
totalheight=$(cat positions.dat | gawk '{print $3}' | LC_ALL=C sort -k 1 -g | tail -n 1)
loc_z=$(echo "$loc_x  $loc_y  $totalheight" | gawk '{print 1.5*sqrt($1*$1+$2*$2+$3*$3)}')
fi

if [ ! -d "Visualization" ]
then
mkdir -p $path/Visualization
fi
rm $path/Visualization/*.dat
# rm $path/Visualization/spinstructure_*.png # make a backup if you don't want to loose all your files!

###########################
# prepare file for povray #
###########################

if [ $s_option == false ] # if no states were given
then
    # states=(start $(seq 1 1 4) end)
    states=(start end)
    # states=(start $(seq 1 1 79) end)
fi
nb_states=${#states[@]}
for state in ${states[@]}
do
if [ $xcut == true ]
then
paste positions.dat magnetic_$state.dat | gawk -v y_cut="$cutpos_y" '{if (strtonum($2) == strtonum(y_cut)) { for (i=1;i<=NF;i++) printf "%.10f%s", $i, (i<NF?OFS:ORS)}}' | gawk '{print "Spin(", "\t", $1",\t",  $2",\t", $3",\t", $4",\t", $5",\t", $6")"}' > Visualization/position_magnetization_$state.dat
# paste positions.dat spin_minimization$state.dat | gawk -v y_cut="$cutpos_y" '{if (strtonum($2) == strtonum(y_cut)) { for (i=1;i<=NF;i++) printf "%.10f%s", $i, (i<NF?OFS:ORS)}}' | gawk '{print "Spin(", "\t", $1",\t",  $2",\t", $3",\t", $4",\t", $5",\t", $6")"}' > Visualization/position_magnetization_$state.dat
else
paste positions.dat magnetic_$state.dat | gawk '{ for (i=1;i<=NF;i++) printf "%.10f%s", $i, (i<NF?OFS:ORS)}' | gawk '{print "Spin(", "\t", $1",\t",  $2",\t", $3",\t", $4",\t", $5",\t", $6")"}' > Visualization/position_magnetization_$state.dat
# paste positions.dat spin_minimization$state.dat | gawk '{ for (i=1;i<=NF;i++) printf "%.10f%s", $i, (i<NF?OFS:ORS)}' | gawk '{print "Spin(", "\t", $1",\t",  $2",\t", $3",\t", $4",\t", $5",\t", $6")"}' > Visualization/position_magnetization_$state.dat
fi
done

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cp $script_dir/povrayscript.pov $path/Visualization
cd $path/Visualization

# sed for linux, gsed for macos
if [ $xcut == true ]
then
gsed -i "/camera {/{n;s/.*/location < $loc_x, $loc_y, $loc_z >\nlook_at  < $loc_x, $cutpos_y, 0 >\nsky < 0, 0, 1 > /}" povrayscript.pov
else
gsed -i "/camera {/{n;s/.*/location < $loc_x, $loc_y, $loc_z >\nlook_at  < $loc_x, $loc_y, 0 > /}" povrayscript.pov
fi
gsed -i "/light_source {/{n;s/.*/        <$loc_x, $loc_y, 25>  \/\/ Location of the source in < x, y, z > /}" povrayscript.pov

if [ $os == true ]
then
    image_num=(${states[@]})
else
    image_num=($(seq -w 0 1 $nb_states))
fi
i=0
for state in ${states[@]}
do
# sed for linux, gsed for macos
gsed -i "$ s/#include .*/#include \"position_magnetization_$state.dat\"/g" povrayscript.pov
##############
# run povray #
##############
echo "generating picture for state: " $state
echo "image number:" ${image_num[$i]}
povray -w${resolution} -h${resolution} povrayscript.pov Output_File_Name=spinstructure_${image_num[$i]}.png &> /dev/null
let i++

done

rm position_magnetization_*.dat

echo "Thank you very much!!!"
echo "Your pictures are now in the directory 'Visualization'"
# you can use the following command in order to generate a movie from these pictures:
# ffmpeg -r 10 -i spinstructure_%2d.png  -c:v libx264 -vf fps=25 -pix_fmt yuv420p anim-.mp4
