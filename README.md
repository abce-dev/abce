<p align="center">
  <img width="500" height="244" src="misc/abce_rectangle.PNG">
</p>

# ABCE: agent-based capacity expansion modeling

ABCE is a module to perform agent-based capacity expansion (CE) modeling for electricity market systems. It simulates decisions made over time by individual profit-motivated entities coexisting in a shared wholesale electricity market environment.

#### Contents
* [Installation](#installation)
  - [Linux / MacOS / Windows Subsystem for Linux](#linux--macos--windows-subsystem-for-linux)
  - [Windows 10](#windows-10)
  - [Optional, Argonne only: installing with A-LEAF](#optional--argonne-only-installing-with-a-leaf)
  - [Optional: installing with CPLEX](#optional-installing-with-cplex)
* [Usage](#usage)
  - [Running ABCE](#running-abce)
  - [Input files](#input-files)
  - [Outputs](#outputs)
* [Contributing](#contributing)
* [Testing](#testing)
* [License](#license)

## Installation

### Linux / MacOS / Windows Subsystem for Linux

1. Clone this repository to your local machine:

    `git clone https://github.com/abce-dev/abce`

2. If using A-LEAF, see the optional [Installing with A-LEAF](#optional--argonne-only-installing-with-a-leaf) section below. Currently, only users with Argonne gitlab credentials can use A-LEAF, but a public release is coming soon!

3. If using CPLEX, see the optional [Installing with CPLEX](#optional-installing-with-cplex) section below

4. Inside your local `abce` directory, run the installation script with:

   `bash ./install.sh`

5. When prompted for the A-LEAF repository, do one of the following:

   * Enter the absolute path to the directory where you cloned A-LEAF, or

   * Press Enter without entering any text, if not using A-LEAF

6. Wait for the installation script to run to completion. Review any errors/issues printed for your reference at the end of execution.

7. Restart your terminal session, or re-source your `.bashrc` file.

8. If using Conda to manage environments, activate the ABCE conda environment with:

   `conda activate abce_env`

9. Rerun the installation script to complete the environment setup:

   `bash ./install.sh`

10. Test the installation using one of the examples:

   `cd examples/single_agent_example`

   If you have CPLEX installed, update the `settings.yml` file in that directory to indicate "CPLEX" as the solver.

   Then run the example case:

   `python ../../run.py --settings_file=./settings.yml --inputs_path=.`

11. Once the previous command runs to completion without failing, generate a precompiled Julia sysimage file from within the `abce/env` directory.

   `julia make_sysimage.jl`


### Windows

1. Download and install [Miniconda](https://docs.conda.io/en/main/miniconda.html)

2. Download and install [Julia 1.8](https://julialang.org/downloads/). Check the box in the installer to add Julia to the PATH.

3. If using A-LEAF, see the optional [Installing with A-LEAF](#optional--argonne-only-installing-with-a-leaf) section below. Currently, only users with Argonne gitlab credentials can use A-LEAF, but a public release is coming soon!

4. If using CPLEX, see the optional [Installing with CPLEX](#optional-installing-with-cplex) section below

5. Using the Anaconda Powershell, clone this repo to your local machine:

   `git clone https://github.com/abce-dev/abce`

6. Create the local conda environment:

   `conda env create -f .\environment_win.yml`

7. Activate the `abce_anv` conda environment:

   `conda activate abce_env`

8. Set the `ABCE_DIR` environment variable to the absolute path to your `abce` repo (e.g. `C:\Users\myname\abce`)

9. Test the installation using one of the examples:

   `cd examples\single_agent_example`

    If you have CPLEX installed, update the `settings.yml` file in that directory to indicate "CPLEX" as the solver.

    Then run the example case:

   `python ..\..\run.py --settings_file=.\settings.yml --inputs_path=.`

10. Once the previous command runs to completion without failing, generate a precompiled Julia sysimage file from within the `abce\env` directory.

   `julia make_sysimage.jl`

### Optional / Argonne only: installing with A-LEAF

1. Clone the ABCE A-LEAF repo:

   `git clone git-out.gss.anl.gov/kbiegel/kb_aleaf`

2. Inside the A-LEAF directory, run the A-LEAF environment setup script:

   `julia make_julia_environment.jl`

3. Test the A-LEAF install by running the following within your A-LEAF directory:

   `julia execute_ALEAF.jl`

4. Once the previous command runs to completion without failing, generate a precompiled Julia sysimage file with:

   `julia make_sysimage.jl`

5. Set the `ALEAF_DIR` environment variable to the absolute path to your A-LEAF repo

### Optional: installing with CPLEX

1. Download the [CPLEX (IBM ILOG STUDIO 20.10) binaries](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing)

2. Run the CPLEX installer, following all instructions.

3. Check that `CPLEX` is installed properly: open the Windows Command Prompt and run the command `$ cplex`. The output should resemble:

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

If the cplex command is not found, try adding the absolute path of your cplex executable to the `$PATH` environment variable with

`$ export CPLEX_STUDIO_BINARIES=/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/`

`$ export PATH=$PATH:$CPLEX_STUDIO_BINARIES`


## Usage

### Running ABCE
ABCE is invoked by the command:

  `python run.py`

run from the top level of the local `abce` directory. This command can accept several optional arguments:

  * `-f`: automatically agree to overwrite any existing database and output files.

  * `--verbosity=k`, where k takes one of the following values:

    * 0: completely silent execution

    * 1: minimal output to mark progression through timesteps and agent turns only

    * 2 (default): basic output showing sub-steps within agent turns

    * 3: maximally verbose, showing results of many calculations and all DEBUG-level messages

  * `--settings_file=<user_file>`: specify the desired settings file with relative or absolute path `<user_file>`. Default: `./settings.yml`

  * `--inputs_path=<user_inputs_dir>`: specify the location of input files, with either a relative or absolute path. Relative path will be relative to the directory from which ABCE is run, not necessarily the directory where the ABCE source code is saved. Default: `./inputs/`

  * `-d`: "demo" mode, pauses execution at the end of each time step to allow the user to review printed outputs

If you'd like to halt execution of ABCE while it is running, press `Ctrl+C` in the terminal session. All data and results generated up to this point will be preserved in the database file (see [Outputs](#outputs)).


### Input files
The input files required to run ABCE are as follows:

 * `settings.yml`: contains all run-specific settings for each simulation. Data specified here supersedes data specified anywhere else.

 * `inputs/`:

   * `agent_specifications.yml`: definitions for the agents: financial parameters, starting portfolios by unit type, and mandatory retirement dates for owned units

   * `C2N_project_definitions.yml`: contains project activity cost and schedule information for coal-to-nuclear projects

   * `demand_data.csv`: normalized peak demand levels per simulated year (used to scale the `peak_demand` parameter)

   * `unit_specs.yml`: construction and operations cost and parameter data for all possible unit types in the model

   * `inputs/ts_data/`:

     * `timeseries_<quantity>_hourly.csv`: hourly timeseries data for each of the following quantities in the system:

       * `load`: normalized to `peak_demand`

       * `wind` and `solar`: wind and solar availability, normalized to the start-of-year installed capacity of each technology, respectively

       * `reg`, `spin`, and `nspin`: ancillary service procurement requirements, in absolute terms (not scaled)

### Outputs
#### `abce/outputs/<scenario_name>/abce_db.db`
The main output from ABCE is the run's database file. This file contains a running ledger of most pieces of data loaded or generated by the simulation as it runs. The database is continually updated as the code runs, so even if the code crashes or you force-quit it with Ctrl+C, the database will retain all data as it existed when execution halted.


#### `abce/outputs/<scenario_name>/outputs.xlsx`
If the code runs successfully to completion, an xlsx dump of the final state of the database is created. No postprocessing is done on the database, but it's easier to open an xlsx file than explore a .db file in most cases, so this is just for convenience.


#### System and agent portfolio evolution charts
ABCE also outputs plots of the evolution of the overall system portfolio and the individual agents' portfolios throughout the simulation. Results are shown for all simulated years, plus three years projected into the future based on already-set construction and retirement plans as of the final simulated year.


### Use ABCE with `watts`
The[ Workflow and Template Toolkit for Simulations (`watts`)](https://github.com/watts-dev/watts) has a plugin for ABCE. Please see the `watts` documentation for usage. This workflow tool is useful for conducting sensitivity analyses and other experiments with ABCE.

## Testing

### Julia Unit Tests

Julia tests may be run with the following command from within the `test/` directory:

`bash run_julia_tests.sh`


## License
Copyright 2023 Argonne National Laboratory

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
