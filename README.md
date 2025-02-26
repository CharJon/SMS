# sms - (Sparse | SCIP | Scientific) MaxCut Solver

## Project

### Solving MaxCut to Optimality

The main goal of this project is to solve the Maximum Cut Problem (MaxCut for short) to proven optimality.

### Authors

Authors: The main author of this project is **Jonas Charfreitag**.

Contributors who worked as student assistants on this project (in alphabetical order):

- [Mohammed Ghannam](https://github.com/mmghannam)
- Claude Jordan
- Michael Kaibel
- [Franziska Kehe](https://github.com/fkehe)
- [Konstantinos Papadopoulos](https://github.com/FreestyleBuild)

## External dependencies

All but one dependency of this project are manged by git (via submodules).
Only SCIP needs to be "installed". The easiest way to make SCIP available on your machine is using one of the "
installers" ZIB offers for download here: https://scipopt.org/index.php#download.

## Building

1) clone project
2) make sure SCIP is installed (see step above)
3) run ```git submodule update --init --recursive``` (see further below for bugfixes)
4) now cmake building should work:

```bash
mkdir build && cd build
cmake ..
make -j 10 # -j 10 makes make use 10 threads
```

If SCIP is not on your path, or you want to use a different version of SCIP as the one on your path:

```bash
mkdir build && cd build
cmake .. -DSCIP_DIR:PATH="{my/path/to/scip}"
make -j 10 # -j 10 makes make use 10 threads
```

If you want to use a specific compiler:

```bash
cmake .. -DCMAKE_C_COMPILER="{/abs/path/to/eg/gcc}" -DCMAKE_CXX_COMPILER="{/abs/path/to/eg/g++}"
```

## Running

The main executable is called `sms`. It can be found in the build directory after building the project.

```bash
./sms --help
```

## Tests

### Run tests via command line

1) build the project as described above (all CMake targets including tests get build by default)
2) stay in the build directory
3) execute the test binary

```bash
./test/test_test --gtest_filter=*
```

`--gtest_filter` can be used to select the test to run. For example `--gtest_filter=IoTest*` will run all test in the
`IoTest` testsuite, ` --gtest_filter=IoTest.ParseEdgelist` will only run the `ParseEdgelist` test.

The `gtest_filter` can be changed to only execute test matching the expression.

### Writing tests

The main CMakeList.txt includes a macro to add tests based on googles gtest library.
The macro is called sms_add_test.
When writing tests, make sure to

- Test relevant special cases, to cover correctness and completeness as good as possible
- Try to avoid as many dependencies as possible, to keep the tests fast and independent

## Code quality

### Style guide

When in doubt, we follow the google style guide for
c++. [Full version](https://google.github.io/styleguide/cppguide.html).
The most important conventions for naming etc. we follow are configured in the ".clang-format" and ".clang-tidy" files!

### Clang tidy

1) build the project as described above, if necessary delete the build directory.
   CMake will automatically create a json file (compile_commands.json) clang-tidy will use in the next step.
3) run clang tidy on the desired files (from inside the project diretory)

```bash
cd ..
clang-tidy --checks=* -p=build src/*
```

The checks parameter can be adjusted, see `clang-tidy --list-checks` all options. You can type `--checks=*` to use all
checks or `--checks=google*` for all google style checks.

The files to check can also be specified (e.g. `src/maxcut.cpp`).

#### Using our own clang-tidy config

Run the following command to use our own clang-tidy config, which is based on google.

```bash
clang-tidy -config-file=.clang-tidy -p=build src/*/*.cpp
```

### clang-format

We also set up a clang-format file. It can be used with CLion as described
here: https://stackoverflow.com/questions/34648255/using-clang-format-in-clion

In the console you can use

```bash
clang-format  -style=file -i src/*/*.cpp include/sms/*/*.hpp
```

to automatically apply the format changes. If you ommit `-i` the changes will be listed but not applied.

## git submodules

git submodules are as awesome as annoying. How to hardreset all
submodules:
```rm -fr .git/config .git/modules && git submodule deinit -f . && git submodule update --init --recursive```

## cli arguments

- GNU suggets: https://www.gnu.org/prep/standards/html_node/Option-Table.html#Option-Table

## SCIP

- SCIP Parameter: https://www.scipopt.org/doc/html/PARAMETERS.php
- SCIP basic
  concepts: http://co-at-work.zib.de/berlin2009/downloads/2009-09-22/2009-09-22-1100-SH-Basic-Concepts-SCIP.pdf
- SCIP 8.0 release notes: https://arxiv.org/pdf/2112.08872.pdf

#### Cutting Planes in SCIP

- Information on how to implement constraint handlers can be found here: https://www.scipopt.org/doc/html/CONS.php
- Useful slides are here: https://www.scipopt.org/download/slides/SCIP-cuttingPlanes.pdf
- Information for the TSP cutting plane example: https://www.scipopt.org/doc/html/TSP_MAIN.php
- CH slides 2020: https://co-at-work.zib.de/slides/Donnerstag_17.9/SCIP.pdf
- CH slides 2015: http://co-at-work.zib.de/files/20150930-Hendel-slides.pdf

### Event Handler in SCIP

- How to add: https://www.scipopt.org/doc/html/EVENT.php
- Example: https://www.scipopt.org/doc/html/EVENTHDLR_MAIN.php
- Events: https://www.scipopt.org/doc/html/type__event_8h.php

### Primal Heuristics in SCIP

- How to add: https://www.scipopt.org/doc/html/HEUR.php

### Visualize B&B-Tree

- https://www.scipopt.org/doc/html/FAQ.php#visualizebranchandbound

#### Memory Management of and in SCIP

- Reference counter for SCIP objects (like variables, constraints, etc.): https://www.scipopt.org/doc/html/OBJ.php
- Allocating memory, the SCIP way: https://www.scipopt.org/doc/html/MEMORY.php

### Compiling SCIP and linking against different LP-solver

Download link: https://scipopt.org/download/release/scipoptsuite-8.0.2.tgz

### Advanced

Not needed by default for this project

#### Ipopt

To allow "-DIPOPT_DIR=$IPOPT_DIR" flag on Linux:

```bash
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build && cd build
../configure --prefix=/path/to/install/dir
make -j
make install
```

#### Compiling SCIP

Download link: https://scipopt.org/download/release/scipoptsuite-9.1.0.tgz

#### With CPLEX

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/opt/scip-X.X.X-cplex -DLPS=cpx -DCPLEX_DIR=$CPLEX_DIR -DIPOPT=false
```
