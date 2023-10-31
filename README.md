# BIOF501A

## Environment Setup
```
# download repository
git clone https://github.com/AmosFong1/BIOF501A

# download KrakenTools
git clone https://github.com/jenniferlu717/KrakenTools

# remove faulty sym link from krona install
rm -rf .conda/envs/BIOF501A/opt/krona/taxonomy

# make directory to store krona database
mkdir -p krona/taxonomy

# create sym link to krona database
ln -s krona/taxonomy .conda/envs/BIOF501A/opt/krona/taxonomy

# download krona database
wget -pO krona/taxonomy/taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# run ktUpdateTaxonomy script
ktUpdateTaxonomy.sh --only-build
```
