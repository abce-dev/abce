#!/bin/bash

#################################################################
# Install script initialization
#################################################################

# Exit this script if a command fails
set -o errexit

# Constants
RC_FILE="$HOME/.bashrc"
CONDA_ENV_FILE="environment.yml"
REQ_FILE="requirements.txt"
JULIA_MAKE_FILE="make_julia_environment.jl"

# Check for command-line arguments
#   -a: pre-specify the ALEAF_DIR absolute path as a command-line argument
#   -n: ignore conda for package management (1=force ignore conda)
while getopts a:n flag
do
    case "${flag}" in
        a) aleaf_dir=${OPTARG};;
        n) no_conda=1;;
    esac
done

#################################################################
# Set values for environment variables
#################################################################

# Set ABCE_DIR to the location of this script
abce_dir=$( dirname -- $( readlink -f -- "$0"; ) )
echo "\$ABCE_DIR will be set to $abce_dir"
export ABCE_DIR=$abce_dir

# Set up ALEAF_DIR
# If specified by a command-line argument, don't prompt the user for the
#   ALEAF_DIR location. The empty string "" is a valid command-line value
#   for this directory.
if [[ ! -n ${aleaf_dir+x} ]]; then
    # If not specified in the command line, request a value from the user
    echo "Please enter the absolute path to the top level of the ALEAF source directory."
    echo "If you do not have ALEAF installed, press Enter to leave this variable empty."
    read aleaf_dir
fi
echo "\$ALEAF_DIR will be set to $aleaf_dir"


#################################################################
# Set up the environment
#################################################################

# Determine whether the script is running in a conda environment
# If conda is installed and available for environment management, use it
if [[ -z "$no_conda" ]] && [[ ! -z $( conda --version | grep -Eo "conda.*[0-9]{1,2}\.[0-9]{1,2}\.[0-9]{1,2}" ) && ! -z $( conda info --envs | grep "\*" ) ]]; then
    echo "conda environment detected; using conda to manage python packages";

    # Check for an appropriate environment spec file, set in $CONDA_ENV_FILE
    #   at the top of this script
    if [[ ! -f "$CONDA_ENV_FILE" ]]; then
        # The conda environment specification (environment.yml) file was not found
        echo "$ABCE_DIR/$CONDA_ENV_FILE not found. Please ensure you have a conda environment specification file in the top level of your ABCE directory.";
        echo "The default environment.yml file is available for download at https://github.com/biegelk/abce.";
        echo "If you do not want to use conda to manage the ABCE environment, rerun this script with the no-conda flag enabled:";
        echo ">$ ./install.sh [other_args] -n 1";
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

    # If julia 1.8 is not the current version of julia, download and install
    #   julia-1.8.2
    if [[ -z $( julia --version | grep "1.8" ) ]]; then
        # If julia --version doesn't return something containing the value "1.8",
        #   then julia 1.8 needs to be installed or set to the default local julia
        if [[ ! -d "$HOME/julia-1.8.2" ]]; then
            if [[ ! -f "$HOME/julia-1.8.2-linux-x86_64.tar.gz" ]]; then
                echo "Downloading Julia 1.8.2...";
                wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.2-linux-x86_64.tar.gz -P "$HOME";
            fi

            echo "Installing Julia 1.8.2..."
            tar zxvf "$HOME/julia-1.8.2-linux-x86_64.tar.gz";
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

# If the .bashrc file doesn't already exist: create it
if [[ ! -f "${RC_FILE}" ]]; then
    echo "No .bashrc file found; creating a new one at ${RC_FILE}";
    touch "${RC_FILE}";
fi

# Create an associative array to allow looping over ABCE environment variables
declare -A env_vars
env_vars=( ["ABCE_DIR"]="$abce_dir" ["ALEAF_DIR"]="$aleaf_dir" )

# Update or append values for each environment variable in the rc file
echo "Updating environment variables in ${RC_FILE}"
for var_name in "${!env_vars[@]}";
do
    # Get a list of all occurrences of $var_name in the bashrc file
    readarray -t var_lines < <(grep --null "${var_name}" "${RC_FILE}")

    if [[ ${#var_lines[@]} == 0 ]]; then
        # If there are no occurrences of this variable in the bashrc file,
        #   simply append a line exporting it to the end of the file
        echo "export ${var_name}=${env_vars[$var_name]}" >> "${RC_FILE}";
    else
        # If this variable does occur at least once in the bashrc file:
        if [[ "${var_lines[( ${#var_lines[@]} - 1 )]}" == "export $var_name=${env_vars[$var_name]}" ]]; then
            # If the last line referencing this variable already has the correct
            #   form, don't update anything
            echo "$var_name is already set to ${env_vars[$var_name]} in ${RC_FILE}.";
        else
            # If the last line referencing this variable sets it to some other
            #   value, append a new export statement with the updated value to 
            #   the end of $RC_FILE
            echo "export ${var_name}=${env_vars[$var_name]}" >> "${RC_FILE}";

            # Alert the user that there may be old settings for these variables
            #   left over in the bashrc file
            echo "Note: ${var_name} already appears in .bashrc. It may be advisable to delete outdated export statements for $var_name.";
        fi
    fi

done

# Ensure that julia-1.8.2 is added to the $PATH such that the `julia` command
#   invokes julia-1.8.2 instead of any other version that may be present
echo "Ensuring Julia 1.8.2 is added to \$PATH in ${RC_FILE}";

# Get a list of all occurrences of 'julia' in the bashrc file
readarray -t julia_lines < <(grep --null -E "export PATH=.*julia.*" "${RC_FILE}")

if [[ ! -z $( echo "${julia_lines[( ${#julia_lines[@]} - 1 )]}" | grep -E "1\.8\.[0-9]{1,2}" ) ]]; then
    # If the last line to add a julia-related value to $PATH already has the
    #   correct form, don't update anything
    echo "julia-1.8.2 is already correctly added to \$PATH in ${RC_FILE}.";
else
    # If .bashrc already has a line adding a different version of julia to the
    #    path, let the user know that this will update which version of Julia
    #    is found globally
    echo "This operation will update the path in ${RC_FILE} where Julia is found globally: the 'julia' command will now invoke julia-1.8.2.";
    echo "If you use Julia on this device for other applications which require a Julia version other than 1.8.2, issues may arise.";

    # Append a new export statement with the updated value to 
    #   the end of $RC_FILE
    echo "export PATH=$HOME/julia-1.8.2/bin/:\$PATH" >> "${RC_FILE}";
fi

#################################################################
# Check for CPLEX
#################################################################
if [[ -z $( which cplex ) ]]; then
    echo "Command 'cplex' not found."
    echo "Either CPLEX is not installed, or you haven't added the location of the 'cplex' binary to the path."
    echo "If you can't install CPLEX, be sure to change the 'solver' setting in settings.yml to a different solver (recommended alternative: 'HiGHS')."
elif [[ -z $( echo $( which cplex) | grep -E "201" ) ]]; then
    echo "ABCE and A-LEAF require CPLEX 20.1, but it appears that you have a different version installed."
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
