conda env config vars set MAMBA_NO_BANNER=1 # suppress oversized banner
mamba create -n ASTRAL_soundings jupyterlab jupytext siphon metview metview-python
mamba install eccodes pdbufr
