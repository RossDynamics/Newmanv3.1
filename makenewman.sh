#!/bin/bash
echo "Performing steps to compile Newman ..."
cd ./src

echo "Compiling newman2D ... "
make -f makefile.2D clean
make -f makefile.2D 

echo "Compiling newman3D ... "
make -f makefile.3D clean
make -f makefile.3D 

echo "Compiling newman4D ... "
make -f makefile.4D clean
make -f makefile.4D 

echo "Finished."
echo "Executables are in ./src"