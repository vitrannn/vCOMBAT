#!/bin/bash
echo -e "changing directory.."
cd "$(dirname "$0")"
echo -e "erasing stuff.."
rm -rf C_Code/bin/*
rm -rf bin/*
rm C_Code/CMakeCache.txt
cd C_Code
echo -e "cmaking.."
cmake .
echo -e "making.."
make
cd ..
echo -e "copying binary.."
cp -a C_Code/bin/* bin/
echo -e "done"