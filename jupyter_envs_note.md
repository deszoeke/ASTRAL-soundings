# A note about environments for running julia and python in Jupyter

Projects with both python and julia notebooks require some care with environments.
Environments and libraries installed in `.julia/conda` by Conda and accessed by PyCall,
which are in turn inside julia environments may not be compatible with 
those in python environments. 

Note python environments modify the paths in the shell, while julia environments are 
specified inside the julia executable. It is possible but inadvisable to run julia from
inside python environments, as these may conflict with the recommended python managed
internally by Conda in `.julia/conda`.

## Jupyter paths

`~/opt/anaconda3/bin` | the first vanilla anaconda jupyter
`~/.julia/conda/3/x86_64/bin/jupyter` | the jupyter for running julia notebooks, esp if they have to use Conda python libraries
`~/mambaforge/envs/ASTRAL_soundings/{bin,etc,shared}/jupyter` | for running python, each project should run in its own mamba/conda environment from the shell

## Mamba environments

In ASTRAL_soundings I use python packages from ECMWF et al. to read
BUFR files. My (shell-level) conda stopped being able to solve the environment, so I switched to mamba.
Mamba is basically the same as conda, only package environment are prechecked and so maybe install faster.
The main thing I needed was to reinstall a new package manager that didn't have a complicated `base` environment, then install
needed packages needed for each project in it's own environment.

Something in these environments didn't get along with .julia/conda. I commented out the mamba/conda initialization
that turns on the base environment in `~/.bash_profile` and put it in `~/.local/bin/mambainit.bsh` and `~/.local/bin/condainit.bsh`.
Another suggestion is to use `mamba deactivate` to get back even out of the base environment.


