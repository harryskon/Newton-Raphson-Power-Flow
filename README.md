# Newton-Raphson-Power-Flow
C++ Newton Raphson power flow code for power systems

This is the first attempt to publish C++ code related to power systems operation.   
The aim is to create a benchmark library of emerging power system algorithms in order to evaluate them in any microprocessor platform. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

What things you need to install the software and how to install them

#### Eigen library - [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

##### Installing

Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.

In order to use [eigen](https://eigen.tuxfamily.org/dox/namespaceEigen.html), you just need to download and extract [eigen](https://eigen.tuxfamily.org/dox/namespaceEigen.html)'s source code. In fact, the header files in the [eigen](https://eigen.tuxfamily.org/dox/namespaceEigen.html) subdirectory are the only files required to compile programs using [eigen](https://eigen.tuxfamily.org/dox/namespaceEigen.html). The header files are the same for all platforms. It is not necessary to use CMake or install anything.

##### Compiling and running your first program

There is no library to link to. The only thing that you need to keep in mind when compiling a program is that the compiler must be able to find the Eigen header files. The directory in which you placed Eigen's source code must be in the include path. With GCC you use the -I option to achieve this, so you can compile the program with a command like this:

```
$ g++ -I /path/to/eigen/ my_program.cpp -o my_program 
```

### Source Download and Compilation

After having all the required dependencies installed, acquire the source code by cloning the git repository:

```
$ git clone https://github.com/harryskon/Newton-Raphson-Power-Flow.git
```
Enter the directory and compile

```
$ make 
```
Important notes:

  * The files are compiled using C++14 standard.
  * Using [godbolt.org](https://gcc.godbolt.org/) it appears that the earilest version to support -std=c++14 is GCC 4.9.0 or Clang 3.5.0. Thus, to use the -std=c++14 flag, update g++/gcc.
  * To update gcc/g++ on Raspberry Pi 2 (wheezy) could be could be to install the g++ 4.9 packages from "Jessie". 

  First bring the current Wheezy up-to-date:
  ```
  $ sudo apt-get update
  $ sudo apt-get upgrade
  ```
  Then edit /etc/apt/sources.list so that you replace the string "wheezy" with "jessie":
  ```
  $ sudo cp /etc/apt/sources.list /etc/apt/sources.list.wheezy
  $ sudo vi /etc/apt/sources.list
  ```
  Now update the package list and install the 4.9 version of GCC/G++:
  ```
  $ sudo apt-get update
  $ sudo apt-get install gcc-4.9 g++-4.9
  ```
  After this revert to the "original" package list:
  ```
  $ sudo cp /etc/apt/sources.list.wheezy /etc/apt/sources.list
  $ sudo apt-get update
  ```
  This leaves the original gcc,g++ in place. Now, to compile with the 4.9 version, then either set the CC and CXX env vars accordingly or invoke the compilers as gcc-4.9 or g++-4.9 explicitly.

  The Makefile of the project has been updated accordingly to use any gcc/g++ version >= 4.9.

## Running the code

The executable after compilation ```nrpowerflow``` takes up to two arguments.

The first argument is the number(#) of IEEE-bus system.  

The second argument (optional) is the number of runs the user would like to run the power flow code. 
If this argument is not given by the user then the code runs only once. 

Example (IEEE 14-bus system):
```
$ ./estimation 14
```

## Authors

* **Harrys Kon** - [Personal Website](https://harrys.fyi/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

