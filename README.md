# Inventory-routing problem (IRP) linear program solver.

A C++ implementation of the IRP linear program using Gurobi's API.

## Prerequisites

* CMake.

* C++17 compiler or an early version.

* GUROBI solver (8 or an early version). Academics can obtain it via this [link](https://www.gurobi.com/downloads/gurobi-optimizer-eula/#Reg "Gurobi's register page").

## Compile and run instructions

Go to the source code folder and to compile type:

```sh
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

for the debug version or simply

```sh
cmake -H. -Bbuild
cmake --build build
```

for the release version.

To run with a configuration file:

```sh
$ ./build/bnp_pmpoc -f [configuration file path]
```

See the "example.cfg" file at the "input" folder for an example of the input configuration file.