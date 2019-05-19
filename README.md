# PairG
Pairing reads on genome graphs

## compile
```sh
git clone --recursive https://github.com/cjain7/PairG.git
cd PairG
mkdir build && cd build

cmake -DCMAKE_CXX_COMPILER=icpc \
-DCMAKE_C_COMPILER=icc \
-DKOKKOS_ARCH=SKX \
-DKOKKOS_ENABLE_OPENMP=TRUE \
-DPROTOBUF_DIR=/work/04042/cjain7/stampede2/Intel/utility/protobuf-3.6.0/local_install ..

make -j4
```
