# Welcome to ONModel2




# Prerequisites

Building ONModel2 requires the following software installed:

* A C++14-compliant compiler
* CMake `>= 3.9`

# Building ONModel2

The following sequence of commands builds ONModel2.
It assumes that your current working directory is the top-level directory
of the freshly cloned repository:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

The build process can be customized with the following CMake variables,
which can be set by adding `-D<var>={ON, OFF}` to the `cmake` call:

# Documentation

ONModel2 *should* provide a documentation.
