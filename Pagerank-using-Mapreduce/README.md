# Pagerank-using-Mapreduce

Two major concepts, both coming from Google: the pagerank algorithm and mapreduce paradigm for problem solving; are implemented here using C++. In this directory, we attempt to parallelise the pagerank algorithm by attempting to solve it as a Mapreduce problem. This was done in the three different ways as below:

## Mapreduce-pagerank using a Mapreduce-C++ library
Implemented mapreduce-pagerank using [mapreduce C++ library](https://github.com/cdmh/mapreduce). The high level program file is named **mr-pr-cpp.cpp**, from which we call the map-reduce library functions. In order to compile the code, following steps are required (Note that **Boost** configurations are involved here):
```sh
$ g++-7 mr-pr-cpp.cpp /usr/lib/x86_64-linux-gnu/libboost_system.a /usr/lib/x86_64-linux-gnu/libboost_iostreams.a /usr/lib/x86_64-linux-gnu/libboost_filesystem.a -pthread -o mr-pr-cpp.o
```

The executable should be runnable as:
```sh
./mr-pr-cpp.o ${filename}.txt -o ${filename}-pr-cpp.txt. 
```

## Mapreduce-pagerank using a custom Mapreduce-MPI library
Implemened a custom mapreduce library with MPI. As all functions in the map-reduce library are not necessary, we implemented only the functions needed for pagerank. Used this library again for mapreduce-pagerank, with high level file named as **mr-pr-mpi.cpp**. To compile, perform the following steps:
```sh
$ mpic++ -o mr-pr-mpi.o mr-pr-mpi.cpp mapreduce.cpp
```

If you look into the directory **custom_library**, then you can find the files of **mapreduce.cpp** and **mapreduce.h**. These files will be required to be in the Cwd in order for the above script to generate binary files. And after that, the executable can be run as:
```sh
$ mpirun -np {total_processors} ./mr-pr-mpi.o {filename}.txt -o {filename}-pr-mpi.txt
```

## Mapreduce-pagerank using an existing Mapreduce-MPI library
Instead of a custom mapreduce library with MPI, we made use of an existing [mapreduce MPI library](https://mapreduce.sandia.gov/). Named the high level file **mr-pr-mpi-base.cpp**. To compile, perform the following steps:
```sh
$ mpic++ -c mr-pr-mpi-base.cpp -o mr-pr-mpi-base.o -I src
$ mpic++ -g -O mr-pr-mpi-base.o src/libmrmpi_mpicc.a -o mr-pr-mpi-base
```

The executable will be runnable as follows:
```sh
$ mpirun -np {total_processors} mr-pr-mpi-base ${filename}.txt -o {filename}-pr-mpi-base.txt
```

## Miscellaneous comments
Note that all the commands above assumed the same directory structure being present. It also assumes that the input file is present in a directory named "test" in the cwd of the code. And there is a directory named "output" where it writes the outputs of the pagerank algorithm execution. The program **mr-pr.cpp** is simply the serial Mapreduce-pagerank implementation.

