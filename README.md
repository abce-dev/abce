[![GitHub Actions build status (Linux)](https://github.com/biegelk/abce/workflows/CI/badge.svg?branch=dev)](https://github.com/biegelk/abce/actions/workflows/python-app.yml)

# abce: agent-based capacity expansion modeling

abce is a module to perform agent-based capacity expansion (CE) modeling for electricity market systems. It is designed to be coupled with the A-LEAF software tool (by Argonne National Laboratory).

#### Contents
* [Installation](#installation)
* [Usage](#usage)
    - [Input Files](#input-files)
    - [Settings File Parameters](#settings-file-parameters)
    - [Agent Settings](#agent-settings)
    - [Using `abce` with `watts`](#use-abce-with-watts)
* [Contributing](#contributing)
* [Testing](#testing)
* [License](#license)

## Installation

Below are roughly the steps required to install and run `abce`. Since `abce` is still in active development, your mileage may vary.

### Unix / Windows Subsystem for Linux

1. Clone this repository to your local machine:

    `git clone https://github.com/biegelk/abce`

2. (optional) Clone the [dedicated A-LEAF repository](https://git-out.gss.anl.gov/kbiegel/kb-aleaf) for `abce`

3. (optional) Download [CPLEX (IBM ILOG STUDIO 20.10) binaries](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing) and install, following the installer's instructions

4. Inside your local `abce` directory, run the installation script with:

   `bash ./install.sh`

5. When prompted for the A-LEAF repository, do one of the following:

   * Enter the absolute path to the directory where you cloned A-LEAF, or

   * Press Enter without entering any text, if not using A-LEAF

6. Wait for the installation script to run to completion. Review any errors/issues printed for your reference at the end of execution.

7. Restart your terminal session, or re-source your `.bashrc` file.

8. If using Conda to manage environments, activate the `abce` conda environment with:

   `conda activate abce_env`

9. Rerun the installation script to complete the environment setup:

   `bash ./install.sh`

### Windows


Windows
1. Download and install [julia](https://julialang.org/downloads/)

2. Download and install [Anaconda](https://www.anaconda.com/products/distribution)

3. Install [Java Virtual Machine](https://linuxhint.com/install-java-ubuntu-22-04/)  with `sudo apt install -y openjdk-18-jre`\*\*

4. Download [CPLEX (IBM ILOG STUDIO 20.10) binaries](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing)

5. Install CPLEX, following instructions in the CPLEX installer

6. Check that `CPLEX` is installed properly: open the Windows Command Prompt and run the command `$ cplex`. The output should resemble:

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

7. Clone the `abce` repository

8. Run `julia make_julia_environment.jl` from your top-level `abce` directory

9. Run `python run.py -f` in the top-level `abce` directory

10. Run `julia make_sysimage.jl`

Optional steps if using `A-LEAF`:

11. Clone [this fork](https://git-out.gss.anl.gov/kbiegel/kb-aleaf) of `A-LEAF`

12. Run `julia make_aleaf_environment.jl` from the `kb-aleaf` directory

13. Run `julia execute_ALEAF.jl`

14. Run `julia make_sysimage.jl`

15. Add an `ALEAF_DIR` environment variable

\* Do **not** use `sudo`. Only `bash ./Anaconda_VERSION_xARCH.sh`

\*\* If you run into an error where it cannot find the java virtual machine execute: `$ which java`, which should output something like `usr/bin/java`.

Then try: `$ sudo bash ./ILOG_COS_20.10_LINUX_X86_64.bin LAX_VM /usr/bin/java`

\*\*\* If the cplex command is not found, try adding the absolute path of your cplex executable to the `$PATH` environment variable with

`$ export CPLEX_STUDIO_BINARIES=/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/`

`$ export PATH=$PATH:$CPLEX_STUDIO_BINARIES`


## Usage

The simplest way to run `abce` is with `python run.py -f`. Where the `-f` flag overwrites pre-existing results. `abce` will ask for permission to overwrite if the flag is left out. This method defaults to a `settings.yml` file in the `abce` directory. However, users may create their own file and specify where the settings are with:

`python run.py -f --settings_file /path/to/your/settings_file.yml`

### Input Files

Here is a description of all the input files that one needs to successfully develop their own `abce` simulation.

### Settings File Parameters

Here is a description of each parameter in the settings file.

### Agent Settings

Here is a description of each parameter in the `gc_params.yml` file that describes agent characteristics.

### Use ABCE with `watts`
The[ Workflow and Template Toolkit for Simulations (`watts`)](https://github.com/watts-dev/watts) has an `abce` plugin. Please see the `watts` documentation for usage. This workflow tool is useful for conducting sensitivity analyses and other experiments with `abce`.

## Contributing

## Testing
`abce` is written in both Julia and Python, there is no unified testing framework at this time. The Python and Julia code may be tested separately.

### Python Unit Tests
Python tests may be run with `pytest` in the top-level directory.

### Julia Unit Tests

**NOTE: Tests are still being added `abce`.**

### Integration Tests

These tests run a complete scenario to verify a working installation.

#### HiGHS Solver
The HiGHS solver is an open-source solver competitive with commercial solvers and much faster than common open-source solvers such as `GLPK` and `Cbc`. To run the `HiGHS` test case:

```bash
$ python run.py -f --settings-file=./test/highs_settings.yml
```
This command should produce the following output.
```bash
>>>(base) C:\Users\samgd\Research\argonne\abce>python run.py -f
Using ATB Year 2020
Existing file at C:\Users\samgd\Research\argonne\abce\solver_test.db deleted.       
Creating a new database file at C:\Users\samgd\Research\argonne\abce\solver_test.db.
Database created in file 'C:\Users\samgd\Research\argonne\abce\solver_test.db'.
using specified value: 3
using specified value: 140
WARNING:root:No match (or multiple matches) found for unit type Wind; setting unit_specs value for Fuel to 0.
WARNING:root:No match (or multiple matches) found for unit type Solar; setting unit_specs value for Fuel to 0.
WARNING:root:No sysimage file found at C:\Users\samgd\Research\argonne\abce\dispatch.so. Execution will proceed, but the dispatch
sub-module may run extremely slowly. If you already have a dispatch sysimage file, please move it to the filename 
{dispatch_sysimage_path}. If you do not have a dispatch sysimage file, please run 'julia make_sysimage.jl --mode=dispatch' in this 
directory.
Start ALEAF scenario reduction algorithm!
==== DONE ! ================================
Agent #202 is taking its turn...
[ Info: Solver is `highs`
Agent #202's turn is complete.

Agent #201 is taking its turn...
[ Info: Solver is `highs`
Agent #201's turn is complete.


All agent turns are complete.

Table of all assets:
Start ALEAF scenario reduction algorithm!
==== DONE ! ================================
Agent #201 is taking its turn...
[ Info: Solver is `highs`
Agent #201's turn is complete.

Agent #202 is taking its turn...
[ Info: Solver is `highs`
Agent #202's turn is complete.


All agent turns are complete.

Table of all assets:
```


## License
Copyright 2022 Argonne National Laboratory

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
