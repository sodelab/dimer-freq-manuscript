
#directories
TOP=`pwd`/../
SINDO=$TOP/tools/sindo

#executables
MKRPT=$SINDO/bin/mkrpt
RPTSX=$SINDO/bin/rpt2sx
SINDO=$SINDO/bin/sindo
EVAL=$TOP/tools/evaluate/eval-pot.x

mkdir -p sindo

if [ 0 == 1 ] ; then
FRAG=monomer
for n in {1..3} ; do

 #navigate to working directory
 mkdir -p sindo/$FRAG
 mkdir -p sindo/$FRAG/${n}MR
 cd sindo/$FRAG/${n}MR

 #copy requisite files
 echo 0 > track
 cp $TOP/step/harm/$FRAG/nwchem.out .
 cp $TOP/step/harm/$FRAG/nmode .
 cp $TOP/step/common/mkrpt_temp.inp mkrpt.write.inp
 cp $TOP/step/harm/$FRAG/sindo.inp .

 #edit files for variable parameters
 sed -i '' 's/$Natoms/3/' mkrpt.write.inp
 sed -i '' 's/$ModeSym/0,0,0,0/' mkrpt.write.inp
 sed -i '' "s/\$MR/${n}/" mkrpt.write.inp
 sed -i '' 's/$Nfree/4/' mkrpt.write.inp
 sed -i '' 's/$Linear/True/' mkrpt.write.inp
 sed -i '' 's/$Nbasis/11,11,11,11/' mkrpt.write.inp

 #write quadrature geometries as inp*
 $MKRPT < mkrpt.write.inp > mkrpt.write.out

 #evaluate inp* geometries in pot directory
 mkdir -p pot; find . -name "inp.*" -exec mv {} pot \; cd pot; 
 > PES
 m=`cat ../track`
 for (( i = 0; i < $m; i++ )); do 
  $EVAL inp.$i >> PES
 done

 #clean up pot directory
 zip inp.zip inp*
 mv PES inp.zip ../. 
 find . -type f -exec rm {} \;
 
 #read PES 
 cd ../; cp mkrpt.write.inp mkrpt.pes.inp
 sed -i '' 's/WRT/PES/' mkrpt.pes.inp
 $MKRPT < mkrpt.pes.inp > mkrpt.out

 #generate pot files
 rm q*.pot
 $RPTSX > rptsx.out

 #generate sindo input
 cat ../../../common/sindo_${FRAG}_tail.inp >> sindo.inp 
 sed -i '' "s/\$MR/${n}/" sindo.inp

 #run sindo
 $SINDO < sindo.inp > sindo.out

 cd ../../../
done

exit

fi



FRAG=dimer-intra
if [ 0 == 1 ] ; then
#navigate to working directory
#mkdir -p sindo/$FRAG
mkdir -p sindo/$FRAG/pot
cd sindo/$FRAG/pot

#copy requisite files
echo 0 > track
cp $TOP/step/harm/$FRAG/nwchem.out .
cp $TOP/step/harm/$FRAG/nmode .
cp $TOP/step/common/mkrpt_temp.inp mkrpt.write.inp

#edit files for variable parameters
sed -i '' 's/$Natoms/6/' mkrpt.write.inp
sed -i '' 's/$ModeSym/0,0,0,0,0,0,0,0/' mkrpt.write.inp
sed -i '' "s/\$MR/3/" mkrpt.write.inp
sed -i '' 's/$Nfree/8/' mkrpt.write.inp
sed -i '' 's/$Linear/False/' mkrpt.write.inp
sed -i '' 's/$Nbasis/11,11,11,11,11,11,11,11/' mkrpt.write.inp

#write quadrature geometries as inp*
$MKRPT < mkrpt.write.inp > mkrpt.write.out

> PES
m=`cat track`
for (( i = 0; i < $m; i++ )); do
 $EVAL inp.$i >> PES
done

#clean up pot directory
find . -name "inp.*" -exec rm {} \;

#read PES 
cp mkrpt.write.inp mkrpt.pes.inp
sed -i '' 's/WRT/PES/' mkrpt.pes.inp
$MKRPT < mkrpt.pes.inp > mkrpt.out

#generate pot files
rm q*.pot
$RPTSX > rptsx.out

fi

cd $TOP/step/sindo/$FRAG/pot

for n in {1..3} ; do

 mkdir $TOP/step/sindo/$FRAG/${n}MR
 cp *pot $TOP/step/sindo/$FRAG/${n}MR/.
 cp $TOP/step/harm/$FRAG/sindo.inp $TOP/step/sindo/$FRAG/${n}MR/. 

 cd ../${n}MR

 #generate sindo input
 cat $TOP/step/common/sindo_${FRAG}_tail.inp >> sindo.inp
 sed -i '' "s/\$MR/${n}/" sindo.inp

 #run sindo
 $SINDO < sindo.inp > sindo.out

 cd ../pot

done

exit
 
 #navigate to working directory
 mkdir -p sindo/$FRAG
 mkdir -p sindo/$FRAG/${n}MR
 cd sindo/$FRAG/${n}MR

 #copy requisite files
 echo 0 > track
 cp $TOP/step/harm/$FRAG/nwchem.out .
 cp $TOP/step/harm/$FRAG/nmode .
 cp $TOP/step/common/mkrpt_temp.inp mkrpt.write.inp
 cp $TOP/step/harm/$FRAG/sindo.inp .

 #edit files for variable parameters
 sed -i '' 's/$Natoms/6/' mkrpt.write.inp
 sed -i '' 's/$ModeSym/0,0,0,0,0,0,0,0/' mkrpt.write.inp
 sed -i '' "s/\$MR/${n}/" mkrpt.write.inp
 sed -i '' 's/$Nfree/8/' mkrpt.write.inp
 sed -i '' 's/$Linear/False/' mkrpt.write.inp
 sed -i '' 's/$Nbasis/11,11,11,11,11,11,11,11/' mkrpt.write.inp

 #write quadrature geometries as inp*
 $MKRPT < mkrpt.write.inp > mkrpt.write.out

 #evaluate inp* geometries in pot directory
 #mkdir -p pot 
 #find . -name "inp.*" -exec mv {} pot/ \;
 
 #cd pot

 > PES
 m=`cat track`
 for (( i = 0; i < $m; i++ )); do 
  $EVAL inp.$i >> PES
 done
 
 #clean up pot directory
 #find . -name "inp.*" -exec zip inp.zip {} \;
 #zip inp.zip inp*
 #mv PES inp.zip ../. 
 #mv PES ../.
 find . -name "inp.*" -exec rm {} \;
 #find . -type f -exec rm {} \;
 
 #read PES 
 #cd ../
 cp mkrpt.write.inp mkrpt.pes.inp
 sed -i '' 's/WRT/PES/' mkrpt.pes.inp
 $MKRPT < mkrpt.pes.inp > mkrpt.out

 #generate pot files
 rm q*.pot
 $RPTSX > rptsx.out

 #generate sindo input
 cat ../../../common/sindo_${FRAG}_tail.inp >> sindo.inp 
 sed -i '' "s/\$MR/${n}/" sindo.inp

 #run sindo
 $SINDO < sindo.inp > sindo.out

 cd ../../../
 exit

exit



