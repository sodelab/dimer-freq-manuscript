

MAVI=../../../../tools/MaVi/Debug/f.mavi

mkdir -p mavi
mkdir -p mavi/monomer

mkdir -p mavi/monomer/1MR
cd mavi/monomer/1MR

cp ../../../harm/monomer/monomer-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/4/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/1/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

mkdir -p mavi/monomer/2MR
cd mavi/monomer/2MR

cp ../../../harm/monomer/monomer-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/4/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/2/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

mkdir -p mavi/monomer/3MR
cd mavi/monomer/3MR

cp ../../../harm/monomer/monomer-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/4/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/3/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

#dimer
mkdir -p mavi/dimer

mkdir -p mavi/dimer/1MR
cd mavi/dimer/1MR

cp ../../../harm/dimer/dimer-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/12/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/1/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
sed -i 's/vci=.t./vci=.f./' mavi.inp
sed -i 's/vscfci=.t./vscfci=.f./' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

mkdir -p mavi/dimer/2MR
cd mavi/dimer/2MR

cp ../../../harm/dimer/dimer-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/12/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/2/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
sed -i 's/vci=.t./vci=.f./' mavi.inp
sed -i 's/vscfci=.t./vscfci=.f./' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

mkdir -p mavi/dimer/3MR
cd mavi/dimer/3MR

cp ../../../harm/dimer/dimer-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/12/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/3/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
sed -i 's/vci=.t./vci=.f./' mavi.inp
sed -i 's/vscfci=.t./vscfci=.f./' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

#dimer-intra
mkdir -p mavi/dimer-intra

mkdir -p mavi/dimer-intra/1MR
cd mavi/dimer-intra/1MR

cp ../../../harm/dimer-intra/dimer-intra-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/8/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/1/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

mkdir -p mavi/dimer-intra/2MR
cd mavi/dimer-intra/2MR

cp ../../../harm/dimer-intra/dimer-intra-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/8/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/2/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../

mkdir -p mavi/dimer-intra/3MR
cd mavi/dimer-intra/3MR

cp ../../../harm/dimer-intra/dimer-intra-001.hs 001.hs
cp ../../../common/mavi_temp.inp mavi.inp
sed -i 's/$Nfree/8/' mavi.inp
sed -i 's/$NBASIS/11/' mavi.inp
sed -i 's/$MR/3/' mavi.inp
sed -i 's/$MAXNRG/20/' mavi.inp
$MAVI < mavi.inp > mavi.out
cd ../../../
