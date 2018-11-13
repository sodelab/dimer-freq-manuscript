
mkdir -p harm
mkdir -p harm/monomer
mkdir -p harm/dimer
mkdir -p harm/dimer-intra
mkdir -p harm/dimer-inter

cp common/mass harm/.
cp common/atomic harm/.

# harmonic monomer
cd harm
../../tools/anharmonic/anharmonic.x -i ../opt/monomer.opt.xyz -l
mv 001.hs monomer-001.hs
cp nmode ../common/nmode-monomer
find . -path './*' -prune -type f -exec mv -f {} monomer \;
cd ../

cp common/mass harm/.
cp common/atomic harm/.

# harmonic dimer
cd harm
../../tools/anharmonic/anharmonic.x -i ../opt/dimer.opt.xyz
mv 001.hs dimer-001.hs
cp nmode ../common/nmode-dimer
find . -path './*' -prune -type f -exec mv -f {} dimer \;
cd ../

cp common/mass harm/.
cp common/atomic harm/.

# intermolecular dimer
cd harm
../../tools/anharmonic/anharmonic.x -i ../opt/dimer.opt.xyz -n 12 -a 1 1 1 1 0 0 0 0 0 0 0 0 
mv 001.hs dimer-inter-001.hs
cp nmode ../common/nmode-dimer-inter
find . -path './*' -prune -type f -exec mv -f {} dimer-inter \;
cd ../

cp common/mass harm/.
cp common/atomic harm/.

# intermolecular dimer
cd harm
../../tools/anharmonic/anharmonic.x -i ../opt/dimer.opt.xyz -n 12 -a 0 0 0 0 1 1 1 1 1 1 1 1
mv 001.hs dimer-intra-001.hs
cp nmode ../common/nmode-dimer-intra
find . -path './*' -prune -type f -exec mv -f {} dimer-intra \;
cd ../

