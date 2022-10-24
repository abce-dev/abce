#!/bin/bash

#################################################################
# Install script initialization
#################################################################

# Exit this script if a command fails
set -o errexit

# Constants
RC_FILE="$HOME/.bashrc"
CONDA_ENV_FILE="./environment.yml"

# Check for command-line arguments
while getopts a: flag
do
    case "${flag}" in
        a) aleaf_dir=${OPTARG};;
    esac
done

#################################################################
# Set up environment variables
#################################################################

# Set ABCE_DIR to the location of this script
ABCE_DIR=$( dirname -- $( readlink -f -- "$0"; ) )

echo "Setting environment variable ABCE_DIR=${ABCE_DIR}"
export ABCE_DIR="${ABCE_DIR}"

# Set up ALEAF_DIR
# If specified by a command-line argument, don't prompt the user for the
#   ALEAF_DIR location
if [[ ! -n ${aleaf_dir+x} ]]; then
    # If not specified in the command line, request a value from the user
    echo "Please enter the absolute path to the top level of the ALEAF source directory."
    echo "If you do not have ALEAF installed, press Enter to leave this variable empty."
    read aleaf_dir
fi

echo "Setting environment variable ALEAF_DIR=$aleaf_dir"
export ALEAF_DIR=$aleaf_dir

# Update the .bashrc file
# If the .bashrc file doesn't already exist: create it
echo "${RC_FILE}"
if [[ ! -f "${RC_FILE}" ]]; then
    echo "No .bashrc file found; creating a new one at ${RC_FILE}";
    touch "${RC_FILE}";
fi

# Create an associative array to allow looping over ABCE environment variables
declare -A env_vars
env_vars=( ["ABCE_DIR"]="$ABCE_DIR" ["ALEAF_DIR"]="$ALEAF_DIR" )

# Update or append values for each environment variable in the rc file
for var_name in "${!env_vars[@]}";
do
    if grep -q "${var_name}" "${RC_FILE}"; then
    # If the form "export $var_name" is found at the start of a line,
    #   replace the rest of the line with the appropriate value
        sed -i "s|^export $var_name=.*|export $var_name=${env_vars[$var_name]}|" "${RC_FILE}";
    else
    # If the form "export $var_name" does not start any lines in the rc file,
    #   append a line at the end of the rc file to export this variable
        echo "export ${var_name}=${env_vars[$var_name]}" >> "${RC_FILE}"
    fi
done

#################################################################
# Set up the Python environment
#################################################################

# Determine whether the script is running in a conda environment
# If conda is installed and available for environment management, use it
if [[ ! -z $( conda --version | grep -Eo "conda.*[0-9][0-9]\.[0-9]\.[0-9]" ) && ! -z $( conda info --envs | grep "\*" ) ]]; then
    echo "conda environment detected; using conda to manage python packages"
    # If a conda environment called abce_env doesn't already exist, create it
    if [[ -z $( conda info --envs | grep "abce_env" ) ]]; then
        echo "No conda environment named abce_env found; creating new environment"
        conda env create -f "${CONDA_ENV_FILE}"
    # If this environment already exists, update it
    else
        echo "Found preexisting conda environment named abce_env; updating it"
        conda env update --prefix ./env --file "${CONDA_ENV_FILE}" --prune
    fi
    conda init bash
    conda activate abce_env
# If conda is not available for environment management, use pip to install
#   packages directly
else
    echo "using pip to manage python packages"
    python3 -m pip install --upgrade pip
    if [[ -f "./requirements.txt" ]]; then
        pip install -r "./requirements.txt";
    fi
fi


#################################################################
# Set up the Julia environment
#################################################################

echo "Setting up the local Julia environment..."
julia make_julia_environment.jl --clean

echo "Julia environment set up."

#################################################################
# Cleanup
#################################################################

# Prevent current value of $aleaf_dir from contaminating future runs of this
#   script in the same terminal session
unset aleaf_dir

echo "All setup completed successfully."
