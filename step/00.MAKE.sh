
#directories
TOP=`pwd`/../
SINDO=$TOP/tools/sindo

cd ../tools

cd routine ; make ; cd ..

cd optimize ; make ; cd ..

cd anharmonic ; make ; cd ..

cd eval ; make ; cd ..

cd mavi/install ; make ; cd ../..


cd sindo ;
sed -i '' "s,\$SINDO_ROOT,${SINDO}," config/make.inc
make ; cd ..

cd nitrogen_v1.9 ; make ; cd ..
