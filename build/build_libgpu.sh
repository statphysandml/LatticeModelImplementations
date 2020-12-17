#!/bin/bash

echo "Not working at the moment"
exit 1

source config.sh

source generate_cmake_lists_txt_file.sh

cd ../
mkdir libgpu
cd libgpu

cmake .. -DCudaUsage="GPU"
make -j9
