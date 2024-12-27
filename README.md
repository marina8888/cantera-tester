Implement a simple volumetric radiation loss. Neglect reabsorption, and assume that all radiation goes into the wall. 

### Getting Started 
This project overwrites the existing Cantera 3.0 with custom radiation models. 
Cantera needs to be installed from source, and recompiled to a new version with each change. 
Instructions for doing this in a conda environment called ct-build.

Includes Cantera as a dependency, howevever in this case cantera is not the main branch code, but rather an internal modified version of the cantera repo. Test files in this repo will not work on the current/main version of cantera and require a reference of the cantera library to point to the internal 'under development' code. 

To compile cantera from source: 

1. create new project folder for testing (calling it an arbitrary prokect name) and inside this folder conda create new environment
2. conda update -n base -c defaults conda
3. conda activate <name of env>
3. git clone <link to github repo>
4. git checkout 3.0 or main
5. install dependencies through anaconda following this url: https://cantera.org/compiling/compilation-reqs.html#sec-compilation-reqs
This part is important! env.yaml install for the file, and pip install everything else.
6. create build using the following config options: 
  scons build -j 8 python_package=full system_yamlcpp=y f90_interface=n skip_slow_tests=yes
7. scons install
8. in test repo create a new seperate folder for testing in the same environment? 

At thie point you will have a project folder and inside that a sub-folder called cantera as a subdirectory. In the project folder, add an empty folder for test scripts as another sub-directory. 
The environment should be in the main folder, but the building (scons build scons install) should take place when you go into the cantera folder only. Also, both folders share the same evironment,which needs to be activated first. Then simply in the test folder write `import cantera as ct`

Please note that should not be installed in the environment using the traditional conda/pip install method because it will interefere with the compiled from source Cantera. 

### Changes and Structure
Current evlauations of radiation and important codes:

<h5>src/OneD/StFlow.cpp/evalResidual</h5> calculates qdotRadiation if m_dot_radiation is True. (Gray body, optically thin model using C2O and H2O only). 
<h5>src/OneD/Sim1D</h5> maps string definitions between new and legacy - not that important. 
<h5>interfaces/cython/onedim.py</h5> radiation_enabled.setter (f.radiation_enabled = True)
<h5>interfaces/cython/_onedim.pxd</h5> enableRadiation and radiationEnabled definitions to pass into C++.
<h5>interfaces/cython/_onedim.pyx</h5> 
<h5> cantera/interfaces/cython/cantera/onedim.py Sensitivity python code