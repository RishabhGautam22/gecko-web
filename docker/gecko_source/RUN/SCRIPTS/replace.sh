#!/bin/sh
# entrées :
# dico : dictionnaire (fort.7 )
# scheme.old : schéma à traiter (fort.17)
# sortie : scheme.new = scheme.old avec les noms raccourcis remplacés par les vrais noms des molécules
echo "Replace code names by formula in scheme"
echo " >>> see scheme.new"
cp ../OUT/dictionary.out dico
cp ../OUT/reactionswithcom.dum scheme.old
sed -i "s/+/ + /g" scheme.old
gawk -f prep1.awk dico scheme.old > scheme.new
rm dico
rm scheme.old

