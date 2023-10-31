# BIOF501A

# Usage
git clone https://github.com/jenniferlu717/KrakenTools
rm -rf ~/.conda/envs/BIOF501A/opt/krona/taxonomy
mkdir -p ~/krona/taxonomy
ln -s ~/krona/taxonomy ~/.conda/envs/BIOF501A/opt/krona/taxonomy
wget -pO ~/krona/taxonomy/taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
ktUpdateTaxonomy.sh --only-build
