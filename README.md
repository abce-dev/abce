# abce: agent-based capacity expansion modeling

abce is a module to perform agent-based capacity expansion (CE) modeling for electricity market systems. It is designed to be coupled with the A-LEAF software tool (by Argonne National Laboratory).

#### Contents
* [Installation](#installation)
* [Usage](#usage)
    - [Input Files](#input-files)
    - [Using `abce` with `watts`](#use-abce-with-watts)
* [Contributing](#contributing)
* [Testing](#testing)
* [License](#license)

## Installation

Below are roughly the steps required to install and run `abce`. Since `abce` is still in active development, your mileage may vary.

|Windows|Unix|
|:----------|:----------|
|Download and Install [julia](https://julialang.org/downloads/)|Download and Install [julia](https://julialang.org/downloads/)|
|Download and Install [Anaconda](https://www.anaconda.com/products/distribution)|Download and Install [Anaconda](https://www.anaconda.com/products/distribution)\*|
||Install [Java Virtual Machine](https://linuxhint.com/install-java-ubuntu-22-04/)  with `sudo apt install -y openjdk-18-jre  `  |
|Download CPLEX (IBM ILOG STUDIO 20.10) Binaries||
|Install CPLEX|Install CPLEX with `sudo bash ./ILOG_COS_20.10_LINUX_X86_64.bin`\*\*|
|Check that `CPLEX` is installed properly by running `$ cplex`\*\*\*|Same as Windows.\*\*\*\*|
|Clone the `abce` repository||
|Run `julia make_julia_environment.jl` from top-level of `abce`||
|Run `python run.py -f`||
|Run `julia make_sysimage.jl`||
|(Optional if using `A-LEAF`)|(Optional if using `A-LEAF`)|
|Install [Gurobi 9.1.2](https://www.gurobi.com/downloads/gurobi-software/)||
|Add the Gurobi executable to `PATH`||
|Clone [this fork](https://git-out.gss.anl.gov/kbiegel/kb-aleaf) of `A-LEAF`||
|Run `julia make_aleaf_environment.jl` from the `kb-aleaf` directory||
|Run `julia execute_ALEAF.jl`||
|Run `julia make_sysimage.jl`||
|Add an `ALEAF_DIR` environment variable||

\* Do **not** use `sudo`. Only `bash ./Anaconda_VERSION_xARCH.sh`

\*\* If you run into an error where it cannot find the java virtual machine execute: `$ which java`, which should output something like `usr/bin/java`.

Then try: `$ sudo bash ./ILOG_COS_20.10_LINUX_X86_64.bin LAX_VM /usr/bin/java`

\*\*\* The output should be
```bash
(base) sdotson@research:~$ cplex

Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 20.1.0.0
  with Simplex, Mixed Integer & Barrier Optimizers
5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
Copyright IBM Corp. 1988, 2020.  All Rights Reserved.

Type 'help' for a list of available commands.
Type 'help' followed by a command name for more
information on commands.
```

\*\*\*\* If the cplex command is not found, try adding the absolute path of your cplex executable to the `$PATH` environment variable with

`$ export CPLEX_STUDIO_BINARIES=/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/`

`$ export PATH=$PATH:$CPLEX_STUDIO_BINARIES`


## Usage

### Input Files


### Use ABCE with `watts`

## Contributing

## Testing
`abce` is written in both Julia and Python, there is no unified testing framework at this time. The Python and Julia code may be tested separately.

### Python Unit Tests
Python tests may be run with
```bash
$ cd abce
$ pytest
```
**NOTE: Tests are still being added `abce`.**


## License
Copyright 2021 Argonne National Laboratory

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
