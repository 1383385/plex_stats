PKG_NAME=plex_stats
USER=tfursten

OS=linux-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION="1.0.0"
conda build ..
conda install plex_stats --use-local