PairG
========================================================================

PairG is a prototype implementation that can be used as a module in a sequence to graph alignment tool. PairG is designed to efficiently validate paired-end distance constraints on sequence graphs (either general or acyclic) in an exact manner. Given a range of allowed distances, and the two vertices, the algorithm computes whether there exists a path of length within the specified range from one vertex to another. The algorithm involves building an index matrix for the complete graph, such that a large set of distance queries can be answered quickly. Existing paired-end read to graph mapping tools have used heuristics, here we have implemented an exact algorithm using sparse linear-algebraic operations. 

## Dependencies
- [cmake](https://cmake.org) version >= 3.1
- A C++ compiler with c++11 support, e.g., GNU `g++` (version 4.9+) or Intel `icpc` (version 15+)
- [Google Protobuf](https://github.com/protocolbuffers/protobuf) library
- Other dependencies (e.g., libvgio, kokkos, catch) are automatically downloaded and compiled 

## Download and compile

The repository and external submodules can be downloaded using the recursive clone.

```sh
git clone --recursive <GITHUB_URL>
```

Next, compile the code using cmake utility:

```sh
mkdir build_directory && cd build_directory
cmake <OPTIONS> ../PairG
make -j4
```

OPTIONS: 
1. `-DPROTOBUF_DIR=<path>` should provide *absolute* path to installation directory of google protobuf library. 
2. `-DKOKKOS_ENABLE_OPENMP=TRUE` option should be specified to enable multi-threading
3. Cmake will automatically look for default C/C++ compilers. To modify the default selection, users can set the two variables `-DCMAKE_CXX_COMPILER=<path to C++ compiler>` and `-DCMAKE_C_COMPILER=<path to C compiler>` if needed. 

After the compilation completes, expect an executable `PairG` in your build\_directory. 

## Usage

* Produce help page
```sh
PairG -h
```

* Build an index, and run 10e6 distance queries between randomly generated vertex pairs. For each query pair, we check if the distance is the range [101,200].
```sh
PairG -m vg -r graph.vg -l 101 -u 200 -c 1000000 -t 24
```

As output, PairG prints the size of input graph, index, and the time it took for indexing and querying.

## Graph input format
PairG accepts a sequence graph (either general or acyclic) in two input formats: `.vg` and `.txt`. `.vg` is a protobuf serialized graph format, defined by VG tool developers [here](https://github.com/vgteam/vg/wiki/File-Formats). `.txt` is a simple human readable format. The first line indicates the count of total vertices (say *n*). Each subsequent line contains information of vertex *i*, 0 <= *i* < *n*. The information in a single line conveys its zero or more out-neighbor vertex ids, followed by its non-empty DNA sequence (either space or tab separated). For example, the following graph is a directed chain of four vertices: `AC (id:0) -> GT (id:1) -> GCCGT (id:2) -> CT (id:3)`

```sh
4
1 AC
2 GT
3 GCCTG
CT
```

The first line above specifies the count of vertices as 4. The second line specifies that vertex 0 has an outgoing edge to vertex 1, and the label of vertex 0 is "AC". Similarly, the third line specifies that vertex 1 has an outgoing edge to vertex 2, and its label is "GT". 
