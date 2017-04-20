#!/bin/bash
echo "Performing steps to compile raw2tec ..."

echo "Compiling tecio library ..."
cd ./lib/tecio
./Runmake
cp tecio.a ../libtecio.a

cd ../../raw2tec
echo "Compiling raw2tec2D ... "
make -f makefile.2D clean
make -f makefile.2D 

echo "Compiling raw2tec3D ... "
make -f makefile.3D clean
make -f makefile.3D 

echo "Compiling raw2tec4D ... "
make -f makefile.4D clean
make -f makefile.4D 

echo "Finished."
echo "Executables are in ./raw2tec."
