#!/bin/bash -i

##########################################################################
# Copyright 2023 Argonne National Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

#################################################################
# Install script initialization
#################################################################

# Exit this script if a command fails
set -o errexit

# Constants
CONDA_ENV_FILE="env/environment_unix.yml"
REQ_FILE="env/requirements.txt"
JULIA_MAKE_FILE="env/make_julia_environment.jl"
JULIA_URL="https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.2-linux-x86_64.tar.gz"

# Check for command-line arguments
#   -n: ignore conda for package management (1=force ignore conda)
#   -f: auto-agree to any yes/no user prompts
while getopts a:nf flag
do
    case "${flag}" in
        n) no_conda=1;;
        f) force=1;;
    esac
done

#################################################################
# Set values for environment variables
#################################################################

# Set ABCE_DIR to the location of this script
abce_dir=$( dirname -- $( readlink -f -- "$0"; ) )
echo "\$ABCE_DIR will be set to $abce_dir"
export ABCE_DIR=$abce_dir
export ABCE_ENV="$ABCE_DIR/env"


#################################################################
# Set up the environment
#################################################################

# Determine whether conda is available; if so, ask the user for permission to 
#   create a dedicated conda environment for ABCE
if [[ -z "$no_conda" ]] && [[ ! -z $( conda --version | grep -Eo "conda.*[0-9]{1,2}\.[0-9]{1,2}\.[0-9]{1,2}" ) && ! -z $( conda info --envs | grep "\*" ) ]] && [[ ! $force ]]; then
    user_resp=""
    while [[ "${user_resp}" != "y" && "${user_resp}" != "n" ]]; do
        echo "I've detected that conda is available on this machine. Can I use conda to create a dedicated environment for abce? (recommended) [y/n]"
        read user_resp
        if [[ $user_resp == "y" ]]; then
            use_conda=true
        fi
    done
fi

# If the user allows use of conda:
if [[ ! -z $use_conda ]]; then
    # Check for an appropriate environment spec file, set in $CONDA_ENV_FILE
    #   at the top of this script
    if [[ ! -f "$CONDA_ENV_FILE" ]]; then
        # The conda environment specification (environment.yml) file was not found
        echo "$ABCE_DIR/$CONDA_ENV_FILE not found. Please ensure you have a conda environment specification file in the top level of your ABCE directory.";
        echo "The default environment_unix.yml file is available for download at https://github.com/biegelk/abce.";
        exit 1;
    else
        # Retrieve the name of the desired conda environment from the yaml file
        CONDA_ENV_NAME=$( grep "name: " "$ABCE_DIR/$CONDA_ENV_FILE" | sed "s|name: ||" );

        # If a conda environment with the name specified in $CONDA_ENV_FILE
        #   doesn't already exist, create it
        if [[ -z $( conda info --envs | grep "$CONDA_ENV_NAME" ) ]]; then
            echo "No conda environment named $CONDA_ENV_NAME found; creating new environment";
            conda env create -f "$ABCE_DIR/$CONDA_ENV_FILE";
        # If this environment already exists, update it
        else
            echo "Found preexisting conda environment named $CONDA_ENV_NAME; updating it"
            conda env update --file "$ABCE_DIR/$CONDA_ENV_FILE" --prune;
        fi
    fi

# If conda is not available for environment management, handle python packages
#   with pip, and other packages via direct download
else
    echo "Using pip to manage python packages"
    python3 -m pip install --upgrade pip;

    # Check for a requirements.txt file in the top-level ABCE directory
    if [[ ! -f "$ABCE_DIR/$REQ_FILE" ]]; then
        echo "$ABCE_DIR/$REQ_FILE not found. Please ensure you have a requirements.txt file in the top level of your ABCE directory.";
        echo "The default requirements.txt file is available for download at https://github.com/biegelk/abce.";
        exit 1;
    else
        pip install -r "$ABCE_DIR/$REQ_FILE";
        echo "All python packages installed with pip.";
    fi

    # Check whether Julia is installed
    if [[ ! -z $( julia --version | grep -E "1\.8\.[0-9]{1,4}" ) ]]; then
        echo "Julia 1.8 is already installed."
    else
        if [[ ! -z $( julia --version | grep "not found" ) ]]; then
            echo "Julia doesn't appear to be installed on this machine."
        elif [[ ! -z $( julia --version | grep -E "[0-9]{1,4}\.[0-9]{1,4}\.[0-9]{1,4}" ) ]]; then
            echo "Your version of Julia is older than 1.8."
        fi

        echo "ABCE requires Julia 1.8."
        user_resp=""
        while [[ "${user_resp}" != "y" || "${user_resp}" != "n" ]]; do
            echo "Can I install Julia 1.8, and set your default Julia version to 1.8? [y/n]"
            read user_resp
        done

        if [[ $user_resp == "y" ]] || [[ $force ]]; then
            # Check for a directory called julia-1.8.2/ in $HOME
            if [[ ! -d "$HOME/julia-1.8.2" ]]; then
                # If there's no such directory, download and unpack the tarball
                #   into $HOME
                if [[ ! -f "$HOME/julia-1.8.2-linux-x86_64.tar.gz" ]]; then
                    # Check whether the julia 1.8.2 tarball has already been
                    #   downloaded into $HOME
                    echo "Downloading Julia 1.8.2..."
                    wget "${JULIA_URL}" -P "$HOME"
                fi

                # Install Julia
                echo "Installing Julia 1.8.2..."
                tar zxvf "$HOME/julia-1.8.2-linux-x86_64.tar.gz"
            fi

        elif [[ $user_resp == "n" ]]; then
            echo "Running ABCE requires Julia 1.8 to be the default, i.e. running 'julia --version' in the local environment returns julia 1.8."
            echo "If you don't want to change the system-wide default, consider installing conda to enable management of multiple package environments."
            echo "You can download Julia 1.8.2 from $JULIA_URL."
            echo "Re-run this installation script once you have Julia 1.8 available as the local default, or once you've installed conda."
            exit 1   
        fi

    fi
fi

# Set up the local Julia environment
echo "Setting up the local Julia environment..."
julia "$ABCE_DIR/$JULIA_MAKE_FILE"
echo "Julia environment created successfully."

echo "ABCE environment created successfully."


#################################################################
# Update the .bashrc file
#################################################################

# Determine the system default shell
shell="${SHELL##*/}"

# Set the configuration file path
RC_FILE="$HOME/.${shell}rc"

# If the shell configuration file doesn't already exist: create it
if [[ ! -f "${RC_FILE}" ]]; then
    echo "No .${shell}rc file found; creating a new one at ${RC_FILE}"
    touch "${RC_FILE}"
fi

# Update or append values for each environment variable in the rc file
echo "Updating environment variables in ${RC_FILE}"

# Create the ABCE update block for bashrc
echo "#==============================================================" >> "${RC_FILE}"
echo "# ABCE configuration" >> "${RC_FILE}"
echo "#   Delete this block to remove undesired side effects (e.g. Julia version update)" >> "${RC_FILE}"
echo "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " >> "${RC_FILE}"
echo "export ABCE_DIR=$abce_dir" >> "${RC_FILE}"
echo "export ABCE_ENV=$abce_dir/env" >> "${RC_FILE}"

# Ensure that julia-1.8.2 is added to the $PATH such that the `julia` command
#   invokes julia-1.8.2 instead of any other version that may be present
echo "Ensuring Julia 1.8.2 is added to \$PATH in ${RC_FILE}"

# Get a list of all occurrences of 'julia' in the bashrc file
IFS=$'\n' julia_lines=($(grep --null -E "export PATH=.*julia.*" "${RC_FILE}" | sed 's/.*=//' | sed 's/".*//'))

if [[ ${#julia_lines[@]} -ne 0 ]]; then
    # If .bashrc already has a line adding a different version of julia to the
    #    path, let the user know that this will update which version of Julia
    #    is found globally
    echo "This operation will update the path in ${RC_FILE} where Julia is found globally: the 'julia' command will now invoke julia-1.8.2."
    echo "If you use Julia on this device for other applications which require a Julia version other than 1.8.2, issues may arise."
    echo "If you didn't want this to happen, open your ${RC_FILE} file and delete the ABCE block at the end of the file."
fi

# Append a new export statement with the updated value to 
#   the end of $RC_FILE
echo "export PATH=$HOME/julia-1.8.2/bin/:\$PATH" >> "${RC_FILE}"

echo "#==============================================================" >> "${RC_FILE}"

#################################################################
# Check for CPLEX
#################################################################
if [[ -z $( which cplex ) ]]; then
    echo "Command 'cplex' not found."
    echo "Either CPLEX is not installed, or you haven't added the location of the 'cplex' binary to the path."
    echo "If you can't install CPLEX, be sure to change the 'solver' setting in settings.yml to a different solver (recommended alternative: 'HiGHS')."
elif [[ -z $( echo $( which cplex) | grep -E "201" ) ]]; then
    echo "ABCE requires CPLEX 20.1, but it appears that you have a different version installed."
    echo "If you can't install CPLEX 20.1, be sure to change the 'solver' setting in settings.yml to a different solver (recommended alternative: 'HiGHS')."
else
    echo "CPLEX 20.1 found on the path."
fi


#################################################################
# Cleanup
#################################################################

# Prevent current values of temporary variables from contaminating future
#   runs of this script in the same terminal session
unset aleaf_dir
unset abce_dir

echo "All setup completed successfully."
echo
echo "Please close and restart your terminal session so the changes can take effect."
if [[ ! -z $CONDA_ENV_NAME ]]; then
    echo "In the new terminal session, switch to the ABCE conda environment with:"
    echo ">$ conda activate $CONDA_ENV_NAME"
fi
echo

