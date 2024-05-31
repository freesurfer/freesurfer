# nnUNet_v1.7 setup instructions

Additional python dependencies are required to run this utility. Follow the
instructions in this file to create a virtual python environment with all the
required dependencies and download the model files required for inference.

## Prerequisites

1. Install conda and locate the path to the activation file, this will be named
'conda.sh'
    * Instructions for installing conda can be found on their site [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
    * Once conda is installed, you can locate the path to your 'conda.sh'
    script by running one of the following commands:
  
    ```bash
    # Linux
    cat ~/.bashrc | grep conda.sh
    ```

    ```zsh
    # MacOS (if conda was initialized for zsh)
    cat ~/.zshrc | grep conda.sh
    ```

2. Download the model file dependencies
    * The model file is available for download [here](https://mega.nz/file/Ph4U0JpA#KMO-hybcDdn_0J5KciI1HO_zYbNPsED_xnHhmCdqoBk)
    * Navigate to the above link and click the 'Download' button in the bottom
    right of the window
    * Once the zip file containing the model is downloaded, unpack it. This can
    be done with the 'Archive Utility' on MacOS, or with the command below on
    either a Linux or MacOS system, after navigating to the directory that
    contains the zip file:

    ```bash
    unzip Task002_StrokesLong.zip
    ```

    * The above command will create unpack the zip file, and create a directory
    named 'Task002_StrokesLong' containing the nnUNet model dependencies. Note
    the location of this directory, as we will need to pass the location to the
    installation script

## Installation script

Once the prerequisites have been met, you can run 'create_nnUNet_v1.7_env.sh' to
setup the nnUNet environment. The script will create a conda env with the tested
version of python and pytorch, clone and install the version of nnUNet used to train the model, create the directory structure to hold the data and model files
for nnUNet, generate a file to source with the env vars required by nnUNet to
find the data.

The following args will need to be passed to the 'create_nnUNet_v1.7_env.sh'
script to correctly configure the nnUNet environment:

* The location of your 'conda.sh' script
* The location of the model file that you have downloaded an unpacked
* The name you want to give the conda env
* The location where you would like to save the trained models and data

Example command:

```bash
./create_env.sh \
    -e linux/etc/profile.d/conda.sh \   # path to conda.sh script
    -m Downloads/Task002_StrokesLong \  # path to unpacked model 
    -n nnUNet_v1.7 \                    # name to give the conda env
    -d nnUNet_paths                     # root of the nnUNet data/model tree
```

When the script finishes running (which could take a bit of time), you will have
a new conda environment with the name passed to the -n flag. A file will also be
created in your pwd, named 'nnUNet_v1.7_path.sh' that contains the env vars that
need to be defined for nnUNet to function properly. Source that file to set the
required env vars.

To have the nnUNet env vars set automatically when the conda env is activated,
move a copy of the nnUNet_v1.7_path.sh into the 'activate.d' directory in your
conda env.

With the newly created conda env active, run the following commands to create
the 'activate.d' directory for the conda env and move a copy of the
'nnUNet_v1.7_path.sh' file there:

```bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
```

```bash
cp nnUNet_v1.7_path.sh $CONDA_PREFIX/etc/conda/activate.d
```
