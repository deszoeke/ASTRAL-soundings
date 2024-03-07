micromamba create -n ASTRAL_soundings jupyterlab jupytext siphon metview metview-python
micromamba activate ASTRAL_soundings
micromamba install -c conda-forge jupyterlab jupytext siphon metview metview-python eccodes pdbufr
