# -- Conda Environmennt --
conda create -n chembl_env python=3.11
conda activate chembl_env
conda install -c conda-forge rdkit
pip install -r requirements.txt

