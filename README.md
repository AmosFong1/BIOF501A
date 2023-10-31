# BIOF501A

## Usage
```
# download repository
git clone https://github.com/AmosFong1/BIOF501A

# create conda environment from environment.yml
conda env create -f environment.yml

# download KrakenTools
git clone https://github.com/jenniferlu717/KrakenTools

# remove faulty sym link from krona install
rm -rf "$(pwd)"/.conda/envs/BIOF501A/opt/krona/taxonomy

# make directory to store krona database
mkdir -p "$(pwd)"/krona/taxonomy

# create sym link to krona database
ln -s "$(pwd)"/krona/taxonomy "$(pwd)"/.conda/envs/BIOF501A/opt/krona/taxonomy

# download krona database
wget -pO "$(pwd)"/krona/taxonomy/taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# run ktUpdateTaxonomy script
ktUpdateTaxonomy.sh --only-build

# run nextflow script
nextflow run BIOF501A.nf
```
