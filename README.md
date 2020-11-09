# Inventory-routing problem (IRP) linear program solver

A C++ implementation of the IRP linear program [[1](#references)] using Gurobi's API.

## Prerequisites

* CMake.

* C++17 compiler or an early version.

* GUROBI solver (9 or an early version). Academics can obtain it via this [link](https://www.gurobi.com/downloads/gurobi-optimizer-eula/#Reg "Gurobi's register page").

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
$ ./build/irp_solver -f [configuration file path]
```

See the "example.cfg" file at the "input" folder for an example of the input configuration file.

## References

**[\[1\] C. Archetti, L. Bertazzi, G. Laporte and M. G. Speranza. A Branch-and-Cut Algorithm for a Vendor-Managed Inventory-Routing Problem Transportation Science, 41(3), 2007, pp. 382-391.](https://pubsonline.informs.org/doi/10.1287/trsc.1060.0188)**