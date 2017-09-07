cd ../../prog/aus/run64/
make -C .. reweight
../reweight
cp pt2gaus.dat ../../../data/pt/pt2gaus_L64.dat
cp pt2gaus.his ../../../data/pt/pt2gaus_L64.his
cp lng.dat     ../../../data/pt/lng_L64.dat
cd ../run64a
../reweight
cp pt2gaus.dat ../../../data/pt/pt2gaus_L64a.dat
cp pt2gaus.his ../../../data/pt/pt2gaus_L64a.his
cp lng.dat     ../../../data/pt/lng_L64a.dat
cd ../run64b
../reweight
cp pt2gaus.dat ../../../data/pt/pt2gaus_L64b.dat
cp pt2gaus.his ../../../data/pt/pt2gaus_L64b.his
cp lng.dat     ../../../data/pt/lng_L64b.dat
cd ../run64c
../reweight
cp pt2gaus.dat ../../../data/pt/pt2gaus_L64c.dat
cp pt2gaus.his ../../../data/pt/pt2gaus_L64c.his
cp lng.dat     ../../../data/pt/lng_L64c.dat
cd ../../../doc/fig
make -B pt_hist.pdf
