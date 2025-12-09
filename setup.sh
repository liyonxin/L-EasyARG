########## Centrifuge
echo ""
echo "Intalling Centrifuge -------------------------------------------------------------------
"
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip --output-document 'centrifuge-1.0.3-beta-Linux_x86_64.zip'
unzip centrifuge-1.0.3-beta-Linux_x86_64.zip
rm centrifuge-1.0.3-beta-Linux_x86_64.zip
rm -rf centrifuge
mv centrifuge-1.0.3-beta ./L-EasyARG-database

########## CARD
echo ""
echo "Intalling CARD -------------------------------------------------------------------
"
wget -c https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2
tar -xjvf broadstreet-v4.0.1.tar.bz2 -C L-EasyARG-database/


##########plasdb
echo ""
echo "Intalling plasdb -------------------------------------------------------------------
"
wget -c https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta
mv download_fasta ./L-EasyARG-database

