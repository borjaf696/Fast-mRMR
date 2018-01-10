# Fast-mRMR
  A shared and distribute memory based fast-mRMR approach. An optimized version for well-known feature selection algorithm mRMR . 
## Compile:

1. ***OpenMP***: 
   * cd Fast-mRMR/fast-mRMR/Código/fast-mRMR-master/cpu-OpenMp/src/Parallel/
   * make clean
   * make
   * export OMP_NUM_THREADS=X
   * ./fast-mrmr -a NumFeatures -f File.mrmr -c ClassNumber
2. ***MPI***  
   * cd Fast-mRMR/fast-mRMR/Código/fast-mRMR-master/cpu-OpenMp/src/Parallel/
   * make clean -f Makefile3
   * make -f Makefile3
   * mpirun -np NumProcs -X OMP_NUM_THREADS=X ./fast-mrmr -a NumFeatures -f File.mrmr -c ClassNumber
   
## mRMR Format:

In order to translate to binary format CSV regular file (comma separed) a conversion tool is offered:
1. ***Data-Reader***: 
  * cd Fast-mRMR/fast-mRMR/Código/fast-mRMR-master/utils/data-reader/
  * make clean
  * make
  * ./mrmr-reader InputFile.csv OutputFile.mrmr

This code is though to be used with CSV which has been previously discretized with at most 256 values and for alphanumeric features with just 1 letter in code. The main reason of that decision is to make it faster enough to deal with large datasets in feasible times.

## Discretization

Data-reader only works properly with datasets which have been previously discretized. A Equal-width based discretizer is also given for python.
  * cd Fast-mRMR/fast-mRMR/Código/fast-mRMR-master/utils/Equal-Width/
  * python Equal.py FileInput FileOutput NumberOfValues NumberOfFeatures Class
  * Char based features will not be discretized because it is supossed to be already bounded to at most 256 different characters.
