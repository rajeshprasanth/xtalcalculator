#!/bin/bash
#
cd ..

if [ ! -d Bin ]; then
mkdir Bin
fi

cd Src 
cp *.x ../Bin/
