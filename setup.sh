conda create -n enz -y
conda init bash
conda activate enz
conda install pip pyrosetta oddt rdkit autodock-vina biopandas scikit-bio -c  https://levinthal:paradox@conda.graylab.jhu.edu -c https://conda.anaconda.org/biocore -c scikit-bio -c oddt -c rdkit -c bioconda -c conda-forge -c anaconda
