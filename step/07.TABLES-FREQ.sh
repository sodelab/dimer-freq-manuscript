#!/bin/bash

#directories
TOP=`pwd`/..

#mkdir supp
#cp common/info.tex supp/info.tex

#functions
round() {
    # adapted from https://stackoverflow.com/questions/26861118/rounding-numbers-with-bc-in-bash
    # $1 is expression to round (should be a valid bc expression)
    # $2 is number of decimal figures (optional). Defaults to two if none given
    if [[ $1 =~ ^[+-]?[0-9]+\.?[0-9]*$ ]]; then 
      local df=${2:-2}
      printf '%.*f\n' "$df" "$(bc -l <<< "a=$1; if(a>0) a+=5/10^($df+1) else if (a<0) a-=5/10^($df+1); scale=$df; a/1")"
    else
      echo $1
    fi
}

##################################################################################################################################################

#
# Extract vibrational frequencies
#

#harmonic frequencies
MONOMER_HARM=( `tail -4 harm/monomer/frequencies.txt | awk '{print $1}'` )
DIMER_HARM=( `tail -12 harm/dimer/frequencies.txt | awk '{print $1}'` )
INTRA_HARM=( `tail -8 harm/dimer-intra/frequencies.txt | awk '{print $1}'` )
INTER_HARM=( `tail -4 harm/dimer-inter/frequencies.txt | awk '{print $1}'` )

#mavi vscf and vmp2 states
MAVI_MONOMER_STATES=( "1000" "0100" "0010" "0001" )
MAVI_MONOMER_VCI_STATES=( "1000" "0100" "0010" "2000" "0001" )
MAVI_INTRA_STATES=( "10000000" "01000000" "00100000" "00010000" "00001000" "00000100" "00000010" "00000001" )
MAVI_INTRA_VCI_STATES=( "10000000" "01000000" "00100000" "00010000" "20000000" "02000000" "00001000" "00000100" "00000010" "00000001" )

#mavi frequencies
MAVI_MONOMER_1MR_VSCF=( `grep -A 5 ' VSCF ' mavi/monomer/1MR/mavi.out | tail -4` )
MAVI_MONOMER_1MR_VMP2=( `grep -A 5 ' VMP2' mavi/monomer/1MR/mavi.out | tail -4` )
MAVI_MONOMER_1MR_VCI=( `head -52 mavi/monomer/1MR/vci.mavi | tail -50 | awk '{print $1}'` )
#MAVI_MONOMER_1MR_VCI=( `grep -A 5 ' VSCFCI' mavi/monomer/1MR/mavi.out | tail -4` )

MAVI_MONOMER_2MR_VSCF=( `grep -A 5 ' VSCF ' mavi/monomer/2MR/mavi.out | tail -4` )
MAVI_MONOMER_2MR_VMP2=( `grep -A 5 ' VMP2' mavi/monomer/2MR/mavi.out | tail -4` )
MAVI_MONOMER_2MR_VCI=( `head -52 mavi/monomer/2MR/vci.mavi | tail -50 | awk '{print $1}'` )
#MAVI_MONOMER_2MR_VCI=( `grep -A 6 ' VSCFCI' mavi/monomer/2MR/mavi.out | tail -5` )

MAVI_MONOMER_3MR_VSCF=( `grep -A 5 ' VSCF ' mavi/monomer/3MR/mavi.out | tail -4` )
MAVI_MONOMER_3MR_VMP2=( `grep -A 5 ' VMP2' mavi/monomer/3MR/mavi.out | tail -4` )
MAVI_MONOMER_3MR_VCI=( `head -52 mavi/monomer/3MR/vci.mavi | tail -50 | awk '{print $1}'` )
#MAVI_MONOMER_3MR_VCI=( `grep -A 6 ' VSCFCI' mavi/monomer/3MR/mavi.out | tail -5` )

#MAVI_DIMER_1MR_VSCF=( `grep -A 13 ' VSCF ' mavi/dimer/1MR/mavi.out` )
#MAVI_DIMER_1MR_VMP2=( `grep -A 13 ' VMP2' mavi/dimer/1MR/mavi.out` )

#MAVI_DIMER_2MR_VSCF=( `grep -A 13 ' VSCF ' mavi/dimer/2MR/mavi.out` )
#MAVI_DIMER_2MR_VMP2=( `grep -A 13 ' VMP2' mavi/dimer/2MR/mavi.out` )

#MAVI_DIMER_3MR_VSCF=( `grep -A 13 ' VSCF ' mavi/dimer/3MR/mavi.out` )
#MAVI_DIMER_3MR_VMP2=( `grep -A 13 ' VMP2' mavi/dimer/3MR/mavi.out` )

MAVI_INTRA_1MR_VSCF=( `grep -A 9 ' VSCF ' mavi/dimer-intra/1MR/mavi.out | tail -8` )
MAVI_INTRA_1MR_VMP2=( `grep -A 9 ' VMP2' mavi/dimer-intra/1MR/mavi.out | tail -8` )
MAVI_INTRA_1MR_VCI=( `head -52 mavi/dimer-intra/1MR/vci.mavi | tail -50 | awk '{print $1}'` )
#MAVI_INTRA_1MR_VCI=( `grep -A 9 ' VSCFCI' mavi/dimer-intra/1MR/vci.mavi | tail -49 | awk '{printf $1}'` )

MAVI_INTRA_2MR_VSCF=( `grep -A 9 ' VSCF ' mavi/dimer-intra/2MR/mavi.out | tail -8 ` )
MAVI_INTRA_2MR_VMP2=( `grep -A 9 ' VMP2' mavi/dimer-intra/2MR/mavi.out | tail -8` )
MAVI_INTRA_2MR_VCI=( `head -52 mavi/dimer-intra/2MR/vci.mavi | tail -50 | awk '{print $1}'` )
#MAVI_INTRA_2MR_VCI=( `grep -A 11 ' VSCFCI' mavi/dimer-intra/2MR/mavi.out | tail -10` )

MAVI_INTRA_3MR_VSCF=( `grep -A 9 ' VSCF ' mavi/dimer-intra/3MR/mavi.out | tail -8` )
MAVI_INTRA_3MR_VMP2=( `grep -A 9 ' VMP2' mavi/dimer-intra/3MR/mavi.out | tail -8` )
MAVI_INTRA_3MR_VCI=( `head -52 mavi/dimer-intra/3MR/vci.mavi | tail -50 | awk '{print $1}'` )
#MAVI_INTRA_3MR_VCI=( `grep -A 11 ' VSCFCI' mavi/dimer-intra/3MR/mavi.out | tail -10` )

#sindo frequencies
SINDO_MONOMER_1MR_VSCF=( `grep "E(VSCF)-E0" sindo/monomer/1MR/sindo.out | awk '{print $NF}'` )
SINDO_MONOMER_1MR_VMP2=( `grep "E(VPT2)-E0" sindo/monomer/1MR/sindo.out | awk '{print $NF}'` )
SINDO_MONOMER_1MR_VCI=( `grep "E(VCI)-E0" sindo/monomer/1MR/sindo.out | awk '{print $NF}'` )

SINDO_MONOMER_2MR_VSCF=( `grep "E(VSCF)-E0" sindo/monomer/2MR/sindo.out | awk '{print $NF}'` )
SINDO_MONOMER_2MR_VMP2=( `grep "E(VPT2)-E0" sindo/monomer/2MR/sindo.out | awk '{print $NF}'` )
SINDO_MONOMER_2MR_VCI=( `grep "E(VCI)-E0" sindo/monomer/2MR/sindo.out | awk '{print $NF}'` )

SINDO_MONOMER_3MR_VSCF=( `grep "E(VSCF)-E0" sindo/monomer/3MR/sindo.out | awk '{print $NF}'` )
SINDO_MONOMER_3MR_VMP2=( `grep "E(VPT2)-E0" sindo/monomer/3MR/sindo.out | awk '{print $NF}'` )
SINDO_MONOMER_3MR_VCI=( `grep "E(VCI)-E0" sindo/monomer/3MR/sindo.out | awk '{print $NF}'` )

SINDO_INTRA_1MR_VSCF=( `grep "E(VSCF)-E0" sindo/dimer-intra/1MR/sindo.out | awk '{print $NF}'` )
SINDO_INTRA_1MR_VMP2=( `grep "E(VPT2)-E0" sindo/dimer-intra/1MR/sindo.out | awk '{print $NF}'` )
SINDO_INTRA_1MR_VCI=( `grep "E(VCI)-E0" sindo/dimer-intra/1MR/sindo.out | awk '{print $NF}'` )

SINDO_INTRA_2MR_VSCF=( `grep "E(VSCF)-E0" sindo/dimer-intra/2MR/sindo.out | awk '{print $NF}'` )
SINDO_INTRA_2MR_VMP2=( `grep "E(VPT2)-E0" sindo/dimer-intra/2MR/sindo.out | awk '{print $NF}'` )
SINDO_INTRA_2MR_VCI=( `grep "E(VCI)-E0" sindo/dimer-intra/2MR/sindo.out | awk '{print $NF}'` )

SINDO_INTRA_3MR_VSCF=( `grep "E(VSCF)-E0" sindo/dimer-intra/3MR/sindo.out | awk '{print $NF}'` )
SINDO_INTRA_3MR_VMP2=( `grep "E(VPT2)-E0" sindo/dimer-intra/3MR/sindo.out | awk '{print $NF}'` )
SINDO_INTRA_3MR_VCI=( `grep "E(VCI)-E0" sindo/dimer-intra/3MR/sindo.out | awk '{print $NF}'` )

#nitrogen frequencies
NITROGEN_INTER_VSCF=( `grep -A 16 '  VSCF  ' nitrogen/vscf/nitrogen-vscf.out | tail -15 | awk '{print $NF}' | tail -14` )
NITROGEN_INTER_VMP2=( `grep -A 16 '  VMP2  ' nitrogen/vmp2/nitrogen-vmp2.out | tail -15 | awk '{print $NF}' | tail -14` )
NITROGEN_INTER_VCI=( `tail -50 nitrogen/vci/nitrogen-vci.out | grep -A 31 'E-E0' | tail -30 | awk '{print $4}' | tail -29` )

##################################################################################################################################################

#
# Extract vibrational states
#

#mavi vscf and vmp2 states
MAVI_MONOMER_1MR_STATES=( "1000" "0100" "0010" "0001" )
MAVI_MONOMER_NMR_STATES=( "1000" "0100" "0010" "0001" )
#MAVI_INTRA_STATES=(10000000 01000000 00100000 00010000 00001000 00000100 00000010 00000001)
#MAVI_INTRA_STATES=(10000000 01000000 00100000 00010000 20000000 02000000 00001000 00000100 00000010 00000001)

#sindo vscf and vmp2 states
SINDO_MONOMER_STATES=( `grep -A 8 '>> VIBRATIONAL STATES' sindo/monomer/1MR/sindo.out | tail -7 | awk '{print $3}'` )
SINDO_INTRA_STATES=( `grep -A 19 '>> VIBRATIONAL STATES' sindo/dimer-intra/1MR/sindo.out | tail -18 | awk '{print $3 $4}'` )

#sindo vci states
SINDO_MONOMER_1MR_STATES=( `grep -A 1000 'ENTER VCI MODULE' sindo/monomer/1MR/sindo.out | grep '> STATE' | awk '{print $4 $5 $6 $7}' | tail -19` )
SINDO_MONOMER_2MR_STATES=( `grep -A 1000 'ENTER VCI MODULE' sindo/monomer/2MR/sindo.out | grep '> STATE' | awk '{print $4 $5 $6 $7}' | tail -19` )
SINDO_MONOMER_3MR_STATES=( `grep -A 1000 'ENTER VCI MODULE' sindo/monomer/3MR/sindo.out | grep '> STATE' | awk '{print $4 $5 $6 $7}' | tail -19` )

MAVI_MONOMER_1MR_STATES=( `head -52 mavi/monomer/1MR/vci.mavi | tail -50 | awk '{print $3 $4 $5 $6}'` )
MAVI_MONOMER_2MR_STATES=( `head -52 mavi/monomer/2MR/vci.mavi | tail -50 | awk '{print $3 $4 $5 $6}'` )
MAVI_MONOMER_3MR_STATES=( `head -52 mavi/monomer/3MR/vci.mavi | tail -50 | awk '{print $3 $4 $5 $6}'` )

MAVI_INTRA_1MR_STATES=( `head -52 mavi/dimer-intra/1MR/vci.mavi | tail -50 | awk '{print $3 $4 $5 $6 $7 $8 $9 $10}'` )
MAVI_INTRA_2MR_STATES=( `head -52 mavi/dimer-intra/2MR/vci.mavi | tail -50 | awk '{print $3 $4 $5 $6 $7 $8 $9 $10}'` )
MAVI_INTRA_3MR_STATES=( `head -52 mavi/dimer-intra/3MR/vci.mavi | tail -50 | awk '{print $3 $4 $5 $6 $7 $8 $9 $10}'` )

SINDO_INTRA_1MR_STATES=( `grep -A 1000 'ENTER VCI MODULE' sindo/dimer-intra/1MR/sindo.out | grep '> STATE' | awk '{print $4 $5 $6 $7 $8 $9 $10 $11}' | tail -49` )
SINDO_INTRA_2MR_STATES=( `grep -A 1000 'ENTER VCI MODULE' sindo/dimer-intra/2MR/sindo.out | grep '> STATE' | awk '{print $4 $5 $6 $7 $8 $9 $10 $11}' | tail -49` )
SINDO_INTRA_3MR_STATES=( `grep -A 1000 'ENTER VCI MODULE' sindo/dimer-intra/3MR/sindo.out | grep '> STATE' | awk '{print $4 $5 $6 $7 $8 $9 $10 $11}' | tail -49` )

NITROGEN_VSCF_STATES=( `grep -A 16 '  VSCF  ' nitrogen/vscf/nitrogen-vscf.out | tail -15 | awk '{print $2 $3 $4 $5}' | tail -14` )
NITROGEN_VMP2_STATES=( `grep -A 16 '  VMP2  ' nitrogen/vmp2/nitrogen-vmp2.out | tail -15 | awk '{print $2 $3 $4 $5}' | tail -14` )
NITROGEN_VCI_STATES=( `tail -50 nitrogen/vci/nitrogen-vci.out | grep -A 31 'E-E0' | tail -30 | awk '{print $7 $8 $9 $10}' | tail -29` )

##################################################################################################################################################


#
# Monomer frequencies
#

declare -a mm1scf; declare -a mm2scf; declare -a mm3scf
declare -a mm1mp2; declare -a mm2mp2; declare -a mm3mp2
declare -a mm1vci; declare -a mm2vci; declare -a mm3vci
declare -a sm1vci; declare -a sm2vci; declare -a sm3vci

for i in "${!SINDO_MONOMER_STATES[@]}"; do 
    mm1scf["${SINDO_MONOMER_STATES[$i]}"]="---"
    mm2scf["${SINDO_MONOMER_STATES[$i]}"]="---"
    mm3scf["${SINDO_MONOMER_STATES[$i]}"]="---"
    mm1mp2["${SINDO_MONOMER_STATES[$i]}"]="---"
    mm2mp2["${SINDO_MONOMER_STATES[$i]}"]="---"
    mm3mp2["${SINDO_MONOMER_STATES[$i]}"]="---"
done

for i in "${!SINDO_MONOMER_3MR_STATES[@]}"; do
    sm1vci["${SINDO_MONOMER_3MR_STATES[$i]}"]="---"
    sm2vci["${SINDO_MONOMER_3MR_STATES[$i]}"]="---"
    sm3vci["${SINDO_MONOMER_3MR_STATES[$i]}"]="---"
    mm1vci["${SINDO_MONOMER_3MR_STATES[$i]}"]="---"
    mm2vci["${SINDO_MONOMER_3MR_STATES[$i]}"]="---"
    mm3vci["${SINDO_MONOMER_3MR_STATES[$i]}"]="---"
done

#for i in "${!MAVI_MONOMER_3MR_STATES[@]}"; do
#    mm1vci["${MAVI_MONOMER_3MR_STATES[$i]}"]="---"
#    mm2vci["${MAVI_MONOMER_3MR_STATES[$i]}"]="---"
#    mm3vci["${MAVI_MONOMER_3MR_STATES[$i]}"]="---"
#done

for i in "${!MAVI_MONOMER_STATES[@]}"; do
    mm1scf["${MAVI_MONOMER_STATES[$i]}"]="${MAVI_MONOMER_1MR_VSCF[$i]}"
    mm2scf[${MAVI_MONOMER_STATES[$i]}]=${MAVI_MONOMER_2MR_VSCF[$i]}
    mm3scf[${MAVI_MONOMER_STATES[$i]}]=${MAVI_MONOMER_3MR_VSCF[$i]}
    mm1mp2[${MAVI_MONOMER_STATES[$i]}]=${MAVI_MONOMER_1MR_VMP2[$i]}
    mm2mp2[${MAVI_MONOMER_STATES[$i]}]=${MAVI_MONOMER_2MR_VMP2[$i]}
    mm3mp2[${MAVI_MONOMER_STATES[$i]}]=${MAVI_MONOMER_3MR_VMP2[$i]}
done

for i in "${!MAVI_MONOMER_3MR_STATES[@]}"; do
    mm1vci[${MAVI_MONOMER_1MR_STATES[$i]}]=${MAVI_MONOMER_1MR_VCI[$i]}
    mm2vci[${MAVI_MONOMER_2MR_STATES[$i]}]=${MAVI_MONOMER_2MR_VCI[$i]}
    mm3vci[${MAVI_MONOMER_3MR_STATES[$i]}]=${MAVI_MONOMER_3MR_VCI[$i]}
done

for i in "${!SINDO_MONOMER_3MR_STATES[@]}"; do
    sm1vci[${SINDO_MONOMER_1MR_STATES[$i]}]=${SINDO_MONOMER_1MR_VCI[$i]}
    sm2vci[${SINDO_MONOMER_2MR_STATES[$i]}]=${SINDO_MONOMER_2MR_VCI[$i]}
    sm3vci[${SINDO_MONOMER_3MR_STATES[$i]}]=${SINDO_MONOMER_3MR_VCI[$i]}
done

mvscf=""
for i in {0..6} ; do
    mvscf+="\$ \\\\nu_{${SINDO_MONOMER_STATES[$i]}} \$ \&"
    mvscf+=" `round ${mm1scf[${SINDO_MONOMER_STATES[$i]}]}` \&"
    mvscf+=" ${SINDO_MONOMER_1MR_VSCF[$i]} \&"
    mvscf+=" `round ${mm2scf[${SINDO_MONOMER_STATES[$i]}]}` \&"
    mvscf+=" ${SINDO_MONOMER_2MR_VSCF[$i]} \&"
    mvscf+=" `round ${mm3scf[${SINDO_MONOMER_STATES[$i]}]}` \&"
    mvscf+=" ${SINDO_MONOMER_3MR_VSCF[$i]} \\\\\\\\ "
done

mvmp2=""
for i in {0..6} ; do
    mvmp2+="\$ \\\\nu_{${SINDO_MONOMER_STATES[$i]}} \$ \&"
    mvmp2+=" `round ${mm1mp2[${SINDO_MONOMER_STATES[$i]}]}` \&"
    mvmp2+=" ${SINDO_MONOMER_1MR_VMP2[$i]} \&"
    mvmp2+=" `round ${mm2mp2[${SINDO_MONOMER_STATES[$i]}]}` \&"
    mvmp2+=" ${SINDO_MONOMER_2MR_VMP2[$i]} \&"
    mvmp2+=" `round ${mm3mp2[${SINDO_MONOMER_STATES[$i]}]}` \&"
    mvmp2+=" ${SINDO_MONOMER_3MR_VMP2[$i]} \\\\\\\\ "
done

#THis is probably OK. But should double check
mvci=""
for i in "${!SINDO_MONOMER_3MR_STATES[@]}"; do
    mvci+="\$ \\\\nu_{${SINDO_MONOMER_3MR_STATES[$i]}} \$ \&"
    mvci+=" `round ${mm1vci[${SINDO_MONOMER_3MR_STATES[$i]}]}` \&"
    mvci+=" `round ${sm1vci[${SINDO_MONOMER_3MR_STATES[$i]}]}` \&"
    mvci+=" `round ${mm2vci[${SINDO_MONOMER_3MR_STATES[$i]}]}` \&"
    mvci+=" `round ${sm2vci[${SINDO_MONOMER_3MR_STATES[$i]}]}` \&"
    mvci+=" `round ${mm3vci[${SINDO_MONOMER_3MR_STATES[$i]}]}` \&"
    mvci+=" `round ${sm3vci[${SINDO_MONOMER_3MR_STATES[$i]}]}` \\\\\\\\ "
done
#BUGBUGBUGBUG

sed -i '' "s,\?monomer-VSCF\?,${mvscf}," supp/info.tex
sed -i '' "s,\?monomer-VMP2\?,${mvmp2}," supp/info.tex
sed -i '' "s,\?monomer-VCI\?,${mvci}," supp/info.tex

##################################################################################################################################################

#
# Intramolecular dimer frequencies
#

declare -a md1scf; declare -a md2scf; declare -a md3scf
declare -a md1mp2; declare -a md2mp2; declare -a md3mp2
declare -a md1vci; declare -a md2vci; declare -a md3vci
declare -a sd1scf; declare -a sd2scf; declare -a sd3scf
declare -a sd1mp2; declare -a sd2mp2; declare -a sd3mp2
declare -a sd1vci; declare -a sd2vci; declare -a sd3vci

for i in "${!SINDO_INTRA_STATES[@]}"; do 
    md1scf["${SINDO_INTRA_STATES[$i]}"]="---"
    md2scf["${SINDO_INTRA_STATES[$i]}"]="---"
    md3scf["${SINDO_INTRA_STATES[$i]}"]="---"
    md1mp2["${SINDO_INTRA_STATES[$i]}"]="---"
    md2mp2["${SINDO_INTRA_STATES[$i]}"]="---"
    md3mp2["${SINDO_INTRA_STATES[$i]}"]="---"
    sd1scf["${SINDO_INTRA_STATES[$i]}"]="---"
    sd2scf["${SINDO_INTRA_STATES[$i]}"]="---"
    sd3scf["${SINDO_INTRA_STATES[$i]}"]="---"
    sd1mp2["${SINDO_INTRA_STATES[$i]}"]="---"
    sd2mp2["${SINDO_INTRA_STATES[$i]}"]="---"
    sd3mp2["${SINDO_INTRA_STATES[$i]}"]="---"
done

for i in "${!SINDO_INTRA_3MR_STATES[@]}"; do 
    md1vci["${SINDO_INTRA_3MR_STATES[$i]}"]="---"
    md2vci["${SINDO_INTRA_3MR_STATES[$i]}"]="---"
    md3vci["${SINDO_INTRA_3MR_STATES[$i]}"]="---"
    sd1vci["${SINDO_INTRA_3MR_STATES[$i]}"]="---"
    sd2vci["${SINDO_INTRA_3MR_STATES[$i]}"]="---"
    sd3vci["${SINDO_INTRA_3MR_STATES[$i]}"]="---"
done

for i in "${!MAVI_INTRA_STATES[@]}"; do
    md1scf[${MAVI_INTRA_STATES[$i]}]=${MAVI_INTRA_1MR_VSCF[$i]}
    md2scf[${MAVI_INTRA_STATES[$i]}]=${MAVI_INTRA_2MR_VSCF[$i]}
    md3scf[${MAVI_INTRA_STATES[$i]}]=${MAVI_INTRA_3MR_VSCF[$i]}
    md1mp2[${MAVI_INTRA_STATES[$i]}]=${MAVI_INTRA_1MR_VMP2[$i]}
    md2mp2[${MAVI_INTRA_STATES[$i]}]=${MAVI_INTRA_2MR_VMP2[$i]}
    md3mp2[${MAVI_INTRA_STATES[$i]}]=${MAVI_INTRA_3MR_VMP2[$i]}
done

for i in "${!SINDO_INTRA_3MR_STATES[@]}"; do
    md1vci[${MAVI_INTRA_1MR_STATES[$i]}]=${MAVI_INTRA_1MR_VCI[$i]}
    md2vci[${MAVI_INTRA_2MR_STATES[$i]}]=${MAVI_INTRA_2MR_VCI[$i]}
    md3vci[${MAVI_INTRA_3MR_STATES[$i]}]=${MAVI_INTRA_3MR_VCI[$i]}
    sd1vci[${SINDO_INTRA_1MR_STATES[$i]}]=${SINDO_INTRA_1MR_VCI[$i]}
    sd2vci[${SINDO_INTRA_2MR_STATES[$i]}]=${SINDO_INTRA_2MR_VCI[$i]}
    sd3vci[${SINDO_INTRA_3MR_STATES[$i]}]=${SINDO_INTRA_3MR_VCI[$i]}
done

for i in "${!MAVI_INTRA_3MR_STATES[@]}"; do
    md1vci[${MAVI_INTRA_1MR_STATES[$i]}]=${MAVI_INTRA_1MR_VCI[$i]}
    md2vci[${MAVI_INTRA_2MR_STATES[$i]}]=${MAVI_INTRA_2MR_VCI[$i]}
    md3vci[${MAVI_INTRA_3MR_STATES[$i]}]=${MAVI_INTRA_3MR_VCI[$i]}
done

dvscf=""
for i in "${!SINDO_INTRA_STATES[@]}"; do
    dvscf+="\$ \\\\nu_{${SINDO_INTRA_STATES[$i]}} \$ \&"
    dvscf+=" `round ${md1scf[${SINDO_INTRA_STATES[$i]}]}` \&"
    dvscf+=" ${SINDO_INTRA_1MR_VSCF[$i]} \&"
    dvscf+=" `round ${md2scf[${SINDO_INTRA_STATES[$i]}]}` \&"
    dvscf+=" ${SINDO_INTRA_2MR_VSCF[$i]} \&"
    dvscf+=" `round ${md3scf[${SINDO_INTRA_STATES[$i]}]}` \&"
    dvscf+=" ${SINDO_INTRA_3MR_VSCF[$i]} \\\\\\\\ "
done

dvmp2=""
for i in "${!SINDO_INTRA_STATES[@]}"; do
    dvmp2+="\$ \\\\nu_{${SINDO_INTRA_STATES[$i]}} \$ \&"
    dvmp2+=" `round ${md1mp2[${SINDO_INTRA_STATES[$i]}]}` \&"
    dvmp2+=" ${SINDO_INTRA_1MR_VMP2[$i]} \&"
    dvmp2+=" `round ${md2mp2[${SINDO_INTRA_STATES[$i]}]}` \&"
    dvmp2+=" ${SINDO_INTRA_2MR_VMP2[$i]} \&"
    dvmp2+=" `round ${md3mp2[${SINDO_INTRA_STATES[$i]}]}` \&"
    dvmp2+=" ${SINDO_INTRA_3MR_VMP2[$i]} \\\\\\\\ "
done

#BUGBUGBUGBUG 
#1MR,2MR,3MR VCI states not necessarily in same order
#THis is probably OK. But should double check
count=0
dvci1=""
dvci2=""
for i in "${!SINDO_INTRA_3MR_STATES[@]}"; do
    if [[ $count -lt 25  ]]; then 
        dvci1+="\$ \\\\nu_{${SINDO_INTRA_3MR_STATES[$i]}} \$ \&"
        dvci1+=" `round ${md1vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci1+=" `round ${sd1vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci1+=" `round ${md2vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci1+=" `round ${sd2vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci1+=" `round ${md3vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci1+=" `round ${sd3vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \\\\\\\\ "
    else
        dvci2+="\$ \\\\nu_{${SINDO_INTRA_3MR_STATES[$i]}} \$ \&"
        dvci2+=" `round ${md1vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci2+=" `round ${sd1vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci2+=" `round ${md2vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci2+=" `round ${sd2vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci2+=" `round ${md3vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \&"
        dvci2+=" `round ${sd3vci[${SINDO_INTRA_3MR_STATES[$i]}]}` \\\\\\\\ "
    fi
    count=$((count+1))
done
#BUGBUGBUGBUG

sed -i '' "s,\?intra-VSCF\?,${dvscf}," supp/info.tex
sed -i '' "s,\?intra-VMP2\?,${dvmp2}," supp/info.tex
sed -i '' "s,\?intra-VCI-1\?,${dvci1}," supp/info.tex
sed -i '' "s,\?intra-VCI-2\?,${dvci2}," supp/info.tex


##################################################################################################################################################

#
# Intermolecular dimer frequencies
#

declare -a nscf; declare -a nmp2; declare -a nvci

for i in "${!NITROGEN_VCI_STATES[@]}"; do 
    nscf["${NITROGEN_VCI_STATES[$i]}"]="---"
    nmp2["${NITROGEN_VCI_STATES[$i]}"]="---"
    nvci["${NITROGEN_VCI_STATES[$i]}"]="---"
done

for i in "${!NITROGEN_VSCF_STATES[@]}"; do
    nscf[${NITROGEN_VSCF_STATES[$i]}]=${NITROGEN_INTER_VSCF[$i]}
done

for i in "${!NITROGEN_VMP2_STATES[@]}"; do
    nmp2[${NITROGEN_VMP2_STATES[$i]}]=${NITROGEN_INTER_VMP2[$i]}
done

for i in "${!NITROGEN_VCI_STATES[@]}"; do
    nvci[${NITROGEN_VCI_STATES[$i]}]=${NITROGEN_INTER_VCI[$i]}
done


count=0
nitro=""
nitro2=""
for i in "${!NITROGEN_VCI_STATES[@]}"; do
    if [[ $count -lt 30  ]]; then 
        nitro+="\$ \\\\nu_{${NITROGEN_VCI_STATES[$i]}} \$ \&"
        nitro+=" `round ${nscf[${NITROGEN_VCI_STATES[$i]}]}` \&"
        nitro+=" `round ${nmp2[${NITROGEN_VCI_STATES[$i]}]}` \&"
        nitro+=" `round ${nvci[${NITROGEN_VCI_STATES[$i]}]}` \\\\\\\\ "
    else
        echo "$count"
        nitro2+="\$ \\\\nu_{${NITROGEN_VCI_STATES[$i]}} \$ \&"
        nitro2+=" `round ${nscf[${NITROGEN_VCI_STATES[$i]}]}` \&"
        nitro2+=" `round ${nmp2[${NITROGEN_VCI_STATES[$i]}]}` \&"
        nitro2+=" `round ${nvci[${NITROGEN_VCI_STATES[$i]}]}` \\\\\\\\  "
    fi
    count=$((count+1))
done

sed -i '' "s,\?inter\?,${nitro}," supp/info.tex
