
#directories
TOP=`pwd`/..
NITROGEN=$TOP/tools/nitrogen_v1.9


#executables
EXE=$NITROGEN/nitrogen/bin/nitrogen
LIB=$TOP/tools/library/libmbCO2-NITROGEN.so

mkdir -p nitrogen
mkdir -p nitrogen/stat
mkdir -p nitrogen/vscf
mkdir -p nitrogen/vmp2
mkdir -p nitrogen/vci

#STAT
cd nitrogen/stat
cp $TOP/step/common/nitrogen-stat.job .
cp $TOP/step/harm/dimer/nitrogen.dat .

ref=`cat nitrogen.dat`

#edit files for variable parameters
sed -i "s,\$LIB,${LIB}," nitrogen-stat.job
sed -i "s,\$REF,${ref}," nitrogen-stat.job

$EXE nitrogen-stat.job > nitrogen-stat.out

cd ../..


#VSCF
cd nitrogen/vscf
cp $TOP/step/common/nitrogen-vib.job nitrogen-vscf.job
cp ../stat/NC_TRANS .

#edit files for variable parameters
sed -i "s,\$LIB,${LIB}," nitrogen-vscf.job

$EXE nitrogen-vscf.job > nitrogen-vscf.out

cd ../..

#VMP2
cd nitrogen/vmp2
cp $TOP/step/common/nitrogen-vib.job nitrogen-vmp2.job
cp ../stat/NC_TRANS .

#edit files for variable parameters
sed -i "s,\$LIB,${LIB}," nitrogen-vmp2.job
echo "# VMP2 parameters" >> nitrogen-vmp2.job
echo "POST_VHF = VMP2"   >> nitrogen-vmp2.job
echo "VCI_MAXV = 10"     >> nitrogen-vmp2.job 

$EXE nitrogen-vmp2.job > nitrogen-vmp2.out

cd ../..


#VCI
cd nitrogen/vci
cp $TOP/step/common/nitrogen-vib.job nitrogen-vci.job
cp ../stat/NC_TRANS .

#edit files for variable parameters
sed -i "s,\$LIB,${LIB}," nitrogen-vci.job
echo "# VCI parameters"  >> nitrogen-vci.job
echo "POST_VHF = VCI"    >> nitrogen-vci.job
echo "VCI_MAXV = 10"     >> nitrogen-vci.job
echo "VCI_NEIG = 30"     >> nitrogen-vci.job

$EXE nitrogen-vci.job > nitrogen-vci.out

cd ../..


