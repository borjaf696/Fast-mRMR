Welcome to the fast-mRMR-MPI wiki!

This is an improved implementation of the classical feature selection method: minimum Redundancy and Maximum Relevance (mRMR); presented by Peng in [1]. 

## Main features

Several optimizations have been introduced in this improved version in order to speed up the costliest computation of the original algorithm: Mutual Information (MI) calculations. These optimizations are described in the followings: 

- **Cache marginal probabilities**: Instead of obtaining the marginal probabilities in each MI computation, those are calculated only once at the beginning of the program, and cached to be re-used in the next iterations.

- **Accumulating redundancy**: The most important optimization is the greedy nature of the algorithm. Instead of computing the mutual information between every pair of features, now redundancy is accumulated in each iteration and the computations are performed between the last selected feature in S and each feature in non-selected set of attributes. 

- **Data access pattern**: The access pattern of mRMR to the dataset is thought to be feature-wise, in contrast to many other ML (machine learning) algorithms, in which access pattern is row-wise. Although being a low-level technical nuance, this aspect can significantly degrade mRMR performance since random access has a much greater cost than block-wise access.

Fast-mRMR-MPI version:
- **Sequential**: original fast-mRMR code for develop features selection.
- **OpenMP**: OpenMP approach allows fast-mRMR to be applied by people or research groups with low resources which leads an acceleration of the original fast-mRMR code scaling with the number of cores.
- **MPI**: MPI approach allows fast-mRMR to be applied over big data datasets in cluster systems.

## Implementations

Here, we include several implementations for different platforms, in order to ease the application of our proposal. These are: 

1. **Sequential version**: we provide a basic implementation in C++ for CPU processing. This is designed to be executed in a single machine. This version includes all aforementioned optimizations.
2. **OpenMP**: an OpenMp's implementation of the algorithm.
3. **MPI**: an Message Passing Interface implementation of the algorithm.


## Project structure

The code is organized as follows:

* _cpu_: C++ code for CPU .
* _cpu-OpenMP: C++ code for shared memory systems.
* _cpu-MPI: C++ code for cluster systems.
* _utils_: this folder contains a data reader program that transforms data in CSV format to the format required by fast-mRMR-MPI algorithm (in binary and columnar-wise format). It also includes a data generator method in case we want to generate synthetic data specifying the structure of this data.
* _examples_: a folder with examples for all versions implemented.   

## Compile:
1. ***Sequential***
   * cd Fast-mRMR/fast-mRMR/cpu/Code/src/
   * make clean
   * make
   * ./fast-mrmr -a NumFeatures -f File.mrmr -c ClassNumber
2. ***OpenMP***: 
   * cd Fast-mRMR/fast-mRMR/cpu-OpenMp/Code/src/
   * make clean
   * make
   * export OMP_NUM_THREADS=X
   * ./fast-mrmr -a NumFeatures -f File.mrmr -c ClassNumber
3. ***MPI***  
   * cd Fast-mRMR/fast-mRMR/cpu-OpenMp/Code/src/
   * make clean
   * make
   * mpirun -np NumProcs -X OMP_NUM_THREADS=X ./fast-mrmr -a NumFeatures -f File.mrmr -c ClassNumber
   
## mRMR Format:

In order to translate to binary format CSV regular file (comma separed) a conversion tool is offered:
1. ***Data-Reader***: 
  * cd Fast-mRMR/fast-mRMR/utils/data-reader/
  * make clean
  * make
  * ./mrmr-reader InputFile.csv OutputFile.mrmr

This code is though to be used with CSV which has been previously discretized with at most 256 values and for alphanumeric features with just 1 letter in code. The main reason of that decision is to make it faster enough to deal with large datasets in feasible times.

## Discretization

Data-reader only works properly with datasets which have been previously discretized. A Equal-width based discretizer is also given for python.
  * cd Fast-mRMR/fast-mRMR/Código/fast-mRMR-master/utils/Equal-Width/
  * python Equal.py FileInput FileOutput NumberOfValues NumberOfFeatures Class
  * Char based features will not be discretized because it is supossed to be already bounded to at most 256 different characters.


