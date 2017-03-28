if [ ! -e path_to_netcdf ]; then
    echo "Path-file path_to_netcdf not found"	
    echo "Please enter directory containing compiled netcdf "
    echo "(usually found in abinit-x.x.x/plugins/netcdf)"
    read -e inputline
    PATH_NETCDF="$inputline"
    if [ ! -e $PATH_NETCDF/netcdf.inc ]; then
	echo "Error: netcdf include files not found at $PATH_NETCDF"
	exit -1
    fi
    echo $PATH_NETCDF > path_to_netcdf
else
    PATH_NETCDF=$(cat path_to_netcdf)       
fi

#Remark for planned transfer of all the codes to abinit-7.x.x: Netcdf seems not to be part of our source codes after abinit version 6

echo
echo
echo Will now compile pspline
echo Including NETCDF from $PATH_NETCDF
echo "(usually found in abinit-x.x.x/plugins/netcdf)"
echo 
echo "proceed (y/n)? "
echo '(if any of the above is incorrect, answer no, delete the appropriate path_to file, and run install again)'
read inputline

what="$inputline"

if [ ! "${what}" = "y" ]; then
    echo Aborting
    exit -1
fi

#Remove old files and build from scratch
#Todo: Introduce the option of simply rebuilding changed files
rm -rf build

# Unpack original code
tar -xzf pspline.tar.gz
# rename folder to "build" to make clear that this is a build folder
# don't make changes here when developping but only in ./patch
mv pspline build

# Create patch
cp ./src_patch/Make.flags ./build/share/
cp ./src_patch/splbrk.f ./build/pspline/




cp $PATH_NETCDF/netcdf.inc build/include
# in the Make.flags file this path is called NETCDF_DIR  
export NETCDF_DIR=$PATH_NETCDF

cd build
make
cd ..
mkdir obj
cd obj
cp ../build/LINUX/lib/*.a .
for l in *.a; do
   ar -x $l
done
ar -r libpsplinetotal.a *.o
cd ..
mv obj/libpsplinetotal.a .
rm -r obj

