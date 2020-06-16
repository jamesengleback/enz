conda create -n enz -y
conda activate enz
conda install -c anaconda pip -y
conda install pyrosetta -c  https://levinthal:paradox@conda.graylab.jhu.edu -y
conda install -c https://conda.anaconda.org/biocore scikit-bio -y
conda install -c oddt oddt -y
conda install rdkit -c rdkit
autodock-vina - docking
conda install -c bioconda autodock-vina
biopandas - cleaning pdbs
conda install -c conda-forge biopandas
