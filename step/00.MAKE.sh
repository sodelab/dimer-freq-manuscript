
#directories
TOP=`pwd`/../
SINDO=$TOP/tools/sindo

cd ../tools

cd routine ; make ; cd ..

cd optimize ; make ; cd ..

cd anharmonic ; make ; cd ..

cd evaluate ; make ; cd ..

#download MaVi from Github
git clone https://github.com/keceli/MaVi.git mavi
cd mavi/Debug ; make ; cd ../..

cd sindo ;
SINDO=$TOP/tools/sindo
sed -i '' "s,\$SINDO_ROOT,${SINDO}," config/make.inc
make ; cd ..

#download NITROGEN from website
cd nitrogen_v1.9 ; make ; cd ..
