# Dev Setup for Samseg Python Port

These instructions apply to a freshly built and upgraded machine running Ubuntu 16.04 and is intended as a reference guide. For a Linux system already set up for develpment, many of these steps can be skipped or appropriately altered.
## Create a dev environment:

```bash
sudo apt-get install build-essential \
            tcsh \
            libtool-bin \
            libtool \
            automake \
            gfortran \
            libglu1-mesa-dev \
            libfreetype6-dev \
            uuid-dev \
            libxmu-dev \
            libxmu-headers \
            libxi-dev \
            libx11-dev \
            libxml2-utils \
            libxt-dev \
            libjpeg62-dev \
            libxaw7-dev \
            liblapack-dev
sudo apt-get install git
```
Check the gcc compiler version:
```bash
gcc --version
```

Be sure to use the same compiler for building both ITK and GEMS2

Most recent _ITK_ requires cmake 3.95 or greater. Test with 

`cmake --version`.

Download and install `cmake-3.10.2-Linux-x86_64.sh`
You may need to copy it to standard location to get it working.
```
sudo cp -r bin /usr/
sudo cp -r doc /usr/share/
sudo cp -r man /usr/share/
sudo cp -r share /usr/
```
## Acquire Source Code
Create a working directory. This document this will refer to it as `~/work/cm`.
```bash
cd ~/work/cm
```
Get the latest release version of ITK:
```bash
git clone https://itk.org/ITK.git
```
Get the porting branch of FreeSurfer and install submodules pybind11 and hdav:
```bash
git clone https://github.com/innolitics/freesurfer.git
cd freesurfer
git checkout  nf-gems2-python-performance
git submodule init
git submodule update
```
## Build from Source
### Build ITK
```bash
cd ~/work/cm/ITK
git checkout release
cd ..
mkdir ITK-build
cd ITK-build
ccmake ../ITK
```
The complete ITK is not needed. You may turn off ```BUILD_TESTING```.

Using the ```t``` command allows you to see the complete set of options.
Then set the ```CMAKE_CXX_FLAGS``` and ```CMAKE_C_FLAGS``` to ```-msse2 -mfpmath=sse```.

The ```c``` command will configure and the ```q``` command will exit.
Now build the ITK library with
```bash
make -j4
```
## Python Setup
Install the python virtual environment tools:
```bash
sudo apt-get install python-pip
sudo apt-get install python3-pip
pip install --upgrade pip
pip install --user virtualenv
pip install --user virtualenvwrapper
```
To take full advantage add this to the end of your `.bashrc` file in your home directory:
```bash
export ITK_DIR=$HOME/work/cm/ITK-build
# where to store our virtual envs
export WORKON_HOME=$HOME/virtenvs
# where projects will reside
export PROJECT_HOME=$HOME/Projects-Active
# where is the virtualenvwrapper.sh
source $HOME/.local/bin/virtualenvwrapper.sh
```
After changing `~/.bashrc` either open a new terminal or do a `source ~/.bashrc`

Locate the python 3 interpreter executable:
```bash
which python3
```
For purposes of this document ```/usr/bin/python3``` will the presumed location.

Now create a `gems2` virtual environment using the correct location of your python3:
```bash
mkvirtualenv gems2 -p /usr/bin/python3
```
The command prompt will indicate you are using this environment with a `(gems2)` prefix.
Having created `gems2` you can enter this environment in any terminal shell by using:
```bash
workon gems2
```
At this point any python packages installed will be local to that environment 
and not interfere with any other python projects, which is the purpose of using virtual environments.

Verify with:
```
which python
```
You should see something like `~/virtenvs/gems2/bin/python3`

### Build and install GEMS2python
In the GEMS2 directory run the ccmake utility
```
cd ~/work/cm/freesurfer/GEMS2
ccmake .
```
Use the `t` option to see the advanced mode options. Then set the `CMAKE_CXX_FLAGS` and `CMAKE_C_FLAGS` to "`-fPIC -fpermissive -msse2 -mfpmath=sse`".

`BUILD_EXECUTABLES` `BUILD_GUI` `BUILD_MATLAB` `BUILD_SHARED_LIBS` and `BUILD_TESTING` should be `OFF`.
`
BUILD_PYTHON` should be `ON`

'CMAKE_BUILD_TYPE' should be 'Release'

Check that the `PYTHON_EXECUTABLE` and `PYTHON_LIBRARY` have valid values such as `/usr/bin/python3.5` 
and `/usr/lib/x86_64-linux-gnu/libpython3.5m.so` respectively.

Build the GEMS2 code with:
```bash
cmake .
make -j4
```
The code just built can be installed into the virtual environment:
```bash
python setup.py bdist_wheel
```
The subdirectory ```dist``` should now have wheel file named something similar to ```GEMS2Python-0.1.0-cp35-cp35m-linux_x86_64.whl```
Install this file into the virtual environment:
```bash
pip install dist/GEMS2Python-0.1.0-cp35-cp35m-linux_x86_64.whl
```

### Build and install samseg
The python portion of the code can also be built and installed as a wheel:
```bash
cd ~/work/cm/freesurfer/samseg/
python setup.py bdist_wheel
pip install dist/samseg-0.1.0-py2.py3-none-any.whl
```
As a test that everything is correctly built, linked, and installed try the following command:
```bash
run_samseg
```
You should see a usage recipe plus an error complaining that no input was specified:
```bash
usage: run_samseg [-h] [-o FOLDER] [-i FILE] [--threads THREADS] [-r FILE]
                  [-m LABEL] [--showfigs] [--nobrainmask] [--diagcovs] [-v]
                  [-a]
run_samseg: error: must specify at least one input
```
## Test and Run

## Running Samseg Code
The matlab script `run_samseg.m` has been ported to `run_samseg_ported.py` with the same command line arguments. 

A minor difference is the use of `-i` and `-o` with single dash instead of `--i` and `--o`
Running with `-h` for help will print a complete usage statement:
```bash
usage: run_samseg_ported.py [-h] [-o FOLDER] [-i FILE] [--threads THREADS]
                            [-r FILE] [--initlta FILE] [-m LABEL] [--movie]
                            [--showfigs] [--nobrainmask] [--diagcovs] [-v]
                            [--reg-only]

optional arguments:
  -h, --help            show this help message and exit
  -o FOLDER, --output FOLDER
                        output to FOLDER
  -i FILE, --input FILE
                        input image(s) from FILE
  --threads THREADS     number of threads
  -r FILE, --regmat FILE
                        skip registration and read from FILE
  --initlta FILE        initial registration FILE
  -m LABEL, --missing LABEL
                        LABEL is a missing structure (repeat for multiple
                        missing labels)
  --movie               show as arrow key controlled time sequence
  --showfigs            show figures during run
  --nobrainmask         no initial brain masking based on affine atlas
                        registration
  --diagcovs            use diagonal covariance matrices (only affect multi-
                        contrast case)
  -v, --verbose         verbose debug output
  --reg-only, --regonly
                        only perform registration
```

### Display Options: showfigs and movie
If you turn on the `showfigs` or `movie` command line options then you will see charted data during the run.

The `showfigs` option by itself will display various graphs and charts as the data is calculated.

The 3d images can be navigated by placing the cursor over one of the three views and then
either clicking with the mouse or using the mouse scroll wheel.
The legend at left also assigns keys to each layer (label or contrast) of the displayed image.
These keys toggle that layer on and off.

The `movie` option alone will show one movie at the end of atas registration and
one movie for each multi resolution level. Time is controlled by using the arrow keys
while the cursor is over the display window.
Left and right arrow keys move forward or back one frame. The up and down arrows move to the
start or end of the movie.
The current frame number and the frame count can be seen in the display title.

If both `showfigs` and `movie` options are selected, then some of the displays will be for a movie in progress,
with as many frames as have been generated to that point in time.

***With either option, all calculations are paused until the displayed window is closed.***

## Running Test Scripts
Place the `innolitics_testing` data folder at `~/work/cm/`


The following will run test case 008 directly:
```bash
cd ~/work/cm/freesurfer
workon gems2
export SAMSEG_DATA_DIR=$HOME/work/cm/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2
run_samseg \
   -o $HOME/work/cm/innolitics_testing/python_temp_data/008 \
    -i $HOME/work/cm/innolitics_testing/buckner40/008/orig.mgz
```

## Development Cycle
If you will be actively developing the code you can uninstall, rebuild, and reinstall the wheels.
But there are some easier alternatives.

### Python Development Cycle
If you will be actively developing the python code then the ```samseg``` wheel can be installed in edit mode.
If you had previously installed it normally, then uninstall it now:
```bash
pip uninstall samseg
```
The project can instead be installed in edit mode.
```bash
cd ~/work/cm/freesurfer/samseg/
pip install -e .
```
Any changes to the python code will be automatically updated.

### C++ Develpment Cycle

A similar process is used if you will also be developing the c++ code.
If you have already installed the ```gems2``` wheel you should uninstall it:
```bash
pip uninstall gems2
```
Instead  add the gems2 binary to the PYTHONPATH:
```bash
export PYTHONPATH="$HOME/work/cm/freesurfer/GEMS2/bin"
```

### The PYTHONPATH alternative
It is also possible to uninstall either or both wheels and use the python path.
```
export PYTHONPATH="$HOME/work/cm/freesurfer/samseg:$HOME/work/cm/freesurfer/GEMS2/bin"
```
Note the use of ```:``` as separator. If you are keeping one of the wheels then remove the corresponding piece of the path.

## Packaging

To create a package to port to another machine change to the main ```samseg``` directory
and run the packaging script ```install_as_virtualenv```
```bash
cd ~/work/cm/freesurfer/samseg/
install_as_virtualenv
```
There are several options. 
The default will create a subfolder ```packaging``` that includes the visualization options.

Running with the ```--help``` option will display a usage message that describes all of the options.
The most significant is the ```-n``` option. This will exclude visualization and
is more suitable for a production server.

The created package is a virtual environment and can be used as such on the build machine.
It includes a copy of the python3 interpreter.
However, the package subfolder can also be copied to another machine with the same OS
and then invoked directly.
It has a script ```run_samseg.sh``` that will use the packaged python
to run the script ```run_samseg``` with the same command line options.
