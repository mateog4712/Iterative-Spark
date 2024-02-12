# Iterative-Spark

#### Description:
Software implementation of Iterative-Spark.      
Iterative-Spark is a time- and space-efficient sparsified minimum free energy folding of possibly-pseudoknotted RNAs

#### Supported OS: 
Linux, macOS


### Installation:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.9.0 or higher), Pthreads, and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that HFold can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/mateog4712/SparseRNAFolD.git) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you need to specify a specific compiler, such as g++, you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   
This can be useful if you are getting errors about your compiler not having C++17 features.

Help
========================================

```
Usage: Iterative-Spark[options] [-r "input structure"] [input sequence]
```

Read input file from cmdline; predict minimum free energy and optimum structure using the time- and space-efficient MFE RNA folding algorithm.

```
  -h, --help             Print help and exit
  -V, --version          Print version and exit
  -v, --verbose          Turn on verbose output
  -m, --mark-candidates  Represent candidate base pairs by curly brackets
  -r, --input-structure  Give a restricted structure as an input structure (required for pseudoknots)
  -p, --pseudoknot-free  Turn off pseudoknot prediction
  -k, --pk-only          Only predict base pairs which cross the given input structure
  -d, --dangles=INT      How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops (default=`2')
  -P, --paramFile        Read energy parameters from paramfile, instead of using the default parameter set.
      --noGC             Turn off garbage collection and related overhead
      --noGU             Turn off Wobble base pairing
```

#### Example:
    assume you are in the Spark directory
    ./build/src/Iterative-Spark GGGGAAAACCCC
    ./build/src/Iterative-Spark -d1 -r "((........))" GGGGAAAACCCC
    ./build/src/Iterative-Spark -d1 -P "src/params/parameters_DP09_Vienna.txt" -r "((........))" GGGGAAAACCCC

```
./build/Iterative-Spark -m -v -r "...............(((((((...)))))))..................." UAACUUAGGGGUUAAAGUUGCAGAUUGUGGCUCUGAAAACACGGGUUCGAA

UAACUUAGGGGUUAAAGUUGCAGAUUGUGGCUCUGAAAACACGGGUUCGAA
{{{({....})}}}.((((((([[[)))))))............]]].... (-7.95)

TA cnt: 63
TA max: 63
TA av:  61
TA rm:  3

Can num:        47
Can cap:        52
TAs num:        63
TAs cap:        64

Psuedoknotted

TAs num:        1
TAs cap:        1
TA av:  5
TA rm:  0

WMB Can num:    6
WMB Can cap:    7
VP Can num:     25
VP Can cap:     25
BE Can num:     7
BE Can cap:     7
```

## Questions
For questions, you can email mateo2@ualberta.ca
