build_dir='build_cmake'

echo "choose compiler option: debug, check, release"
echo " nothing chooses release"
read compop
if [ ${#compop} -eq "0" ]
then compop="release"
fi
echo "using compile option COMPOP: $compop"

githash=`git log --pretty=format:'%H' -n 1`
gitversion=` echo "CPP_VERSIONGIT=\"'$githash'\" "`
if [ -d $build_dir ]; then
    rm -r $build_dir
fi
mkdir -p $build_dir
cp config.cmake $build_dir
cd $build_dir
echo -e "\n\n\n start cmake"
cmake .. -DGITVERSION:STRING=$gitversion -DCOMPOP:STRING=$compop
echo -e "\n\n\n start make"

if command -v nproc &> /dev/null
then
	make VERBOSE=1 -j $(( `nproc` < 8 ? `nproc` : 8 ))
else
	make VERBOSE=1 
fi

make
echo "compop $compop" >> compop
echo -e "\n\n\n done"
bindir=`pwd`
cd ..
ln -s ${bindir}/Matjes . 2> /dev/null
