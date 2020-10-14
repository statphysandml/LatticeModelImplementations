#!/bin/bash

source config.sh

source generate_cmake_lists_txt_file.sh

cd ../
mkdir libgpu
cd libgpu

cmake .. -DCudaUsage="GPU"
make -j9
