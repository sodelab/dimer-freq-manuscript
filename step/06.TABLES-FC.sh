#!/bin/bash

#directories
TOP=`pwd`/..

#mkdir supp
cp common/info.tex supp/info.tex
#cp common/info.tex info.tex

#functions
round() {
    # from https://unix.stackexchange.com/questions/265119/how-to-format-floating-point-number-with-exactly-2-significant-digits-in-bash
    # $1 is precision (significant figures)
    # $2 is expression to round
    n=$(printf "%.${1}g" "$2")
    if [ "$n" != "${n#*e}" ]
    then
        f="${n##*e-}"
        test "$n" = "$f" && f= || f=$(( ${f#0}+$1-1 ))
        printf "%0.${f}f" "$n"
    else
        printf "%s" "$n"
    fi
}

THRESHOLD=0.001

#
# Extract from 001.hs files
#

mfile=mavi/monomer/3MR/001.hs
dfile=mavi/dimer-intra/3MR/001.hs

#1MR monomer Gradient
M1G=( `grep -A 50 "# 1MR" $mfile | grep -A 4 "# Gradient / hartree" | tail -4 | awk '{print $2}'` )
#1MR monomer Hessian
M1H=( `grep -A 50 "# 1MR" $mfile | grep -A 4 "Hessian(i,i)" | tail -4 | awk '{print $2}'` )
#1MR monomer Cubic
M1C=( `grep -A 50 "# 1MR" $mfile | grep -A 4 "Cubic(i,i,i)" | tail -4 | awk '{print $2}'` )
#1MR monomer Quadratic
M1Q=( `grep -A 50 "# 1MR" $mfile | grep -A 4 "Quartic(i,i,i,i)" | tail -4 | awk '{print $2}'` )

#2MR monomer Hessian
M2H=( `grep -A 180 "# 2MR" $mfile | grep -A 6 "Hessian(i,j)" | tail -6 | awk '{print $3}'` )
#2MR monomer Quartic
M2Q1=( `grep -A 180 "# 2MR" $mfile | grep -A 6 "Quartic(i,i,j,j)" | tail -6 | awk '{print $3}'` )
#2MR monomer Cubic
M2C=( `grep -A 180 "# 2MR" $mfile | grep -A 12 "Cubic(i,i,j)" | tail -12 | awk '{print $3}'` )
#2MR monomer Quartic
M2Q2=( `grep -A 180 "# 2MR" $mfile | grep -A 12 "Quartic(i,i,i,j)" | tail -12 | awk '{print $3}'` )

#3MR monomer Cubic
M3C=( `grep -A 250 "# 3MR" $mfile | grep -A 4 "Cubic(i,j,k)" | tail -4 | awk '{print $4}'` )
#3MR monomer Quartic
M3Q=( `grep -A 250 "# 3MR" $mfile | grep -A 12 "Quartic(i,i,j,k)" | tail -12 | awk '{print $4}'` )


#1MR dimer Gradient
D1G=( `grep -A 50 "# 1MR" $dfile | grep -A 8 "# Gradient / hartree" | tail -8 | awk '{print $2}'` )
#1MR dimer Hessian
D1H=( `grep -A 50 "# 1MR" $dfile | grep -A 8 "Hessian(i,i)" | tail -8 | awk '{print $2}'` )
#1MR dimer Cubic
D1C=( `grep -A 50 "# 1MR" $dfile | grep -A 8 "Cubic(i,i,i)" | tail -8 | awk '{print $2}'` )
#1MR dimer Quadratic
D1Q=( `grep -A 50 "# 1MR" $dfile | grep -A 8 "Quartic(i,i,i,i)" | tail -8 | awk '{print $2}'` )

#2MR dimer Hessian
D2H=( `grep -A 180 "# 2MR" $dfile | grep -A 28 "Hessian(i,j)" | tail -28 | awk '{print $3}'` )
#2MR dimer Quartic
D2Q1=( `grep -A 180 "# 2MR" $dfile | grep -A 28 "Quartic(i,i,j,j)" | tail -28 | awk '{print $3}'` )
#2MR dimer Cubic
D2C=( `grep -A 180 "# 2MR" $dfile | grep -A 56 "Cubic(i,i,j)" | tail -56 | awk '{print $3}'` )
#2MR dimer Quartic
D2Q2=( `grep -A 180 "# 2MR" $dfile | grep -A 56 "Quartic(i,i,i,j)" | tail -56 | awk '{print $3}'` )

#3MR dimer Cubic
D3C=( `grep -A 250 "# 3MR" $dfile | grep -A 56 "Cubic(i,j,k)" | tail -56 | awk '{print $4}'` )
#3MR dimer Quartic
D3Q=( `grep -A 250 "# 3MR" $dfile | grep -A 168 "Quartic(i,i,j,k)" | tail -168 | awk '{print $4}'` )


#
# Load into tables
#

mon=""
col=0
for i in {0..3} ; do
    n=$((i+1))

    x=${M1H[$i]}
    if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
        mon+="\$ h_{$n$n} \$ \&"
        mon+=" \$ `round 5 ${M1H[$i]}` \$ \&"
        mon+=" \$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-2}\\\\ {\\\\rm u}^{-1} \$ "
        if [ $((col%2)) -eq 0 ] ; then 
            mon+=" \& ";
        else
            mon+=" \\\\\\\\ ";
        fi
        col=$((col+1))
    fi
done

for i in {0..3} ; do
    n=$((i+1))

    x=${M1C[$i]}
    if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
        mon+="\$ t_{$n$n$n} \$ \&"
        mon+=" \$ `round 5 ${M1C[$i]}` \$ \&"
        mon+=" \$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
        if [ $((col%2)) -eq 0 ] ; then 
            mon+=" \& ";
        else
            mon+=" \\\\\\\\ ";
        fi
        col=$((col+1))
    fi
done

for i in {0..3} ; do
    n=$((i+1))

    x=${M1Q[$i]}
    if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
        mon+="\$ u_{$n$n$n$n} \$ \&"
        mon+=" \$ `round 5 ${M1Q[$i]}` \$ \&"
        mon+=" \$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
        if [ $((col%2)) -eq 0 ] ; then 
            mon+=" \& ";
        else
            mon+=" \\\\\\\\ ";
        fi
        col=$((col+1))
    fi
done

#2MR
count=0
for (( i=0; i<=3; i++ )) ; do
    m=$((i+1))
    for (( j=$m; j<=3; j++ )) ; do
        n=$((j+1))
        
        x=${M2H[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            mon+="\$ h_{$n$m} \$ \&"
            mon+=" \$ `round 5 ${M2H[$count]}` \$ \&"
            mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-2}\\\\ {\\\\rm u}^{-1} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                mon+=" \& ";
            else
                mon+=" \\\\\\\\ ";
            fi
            col=0
        fi
        count=$((count+1))

    done
done

count=0
for (( i=0; i<=3; i++ )) ; do
    m=$((i+1))
    for (( j=$m; j<=3; j++ )) ; do
        n=$((j+1))

        x=${M2Q1[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            mon+="\$ u_{$n$n$m$m} \$ \&"
            mon+=" \$ `round 5 ${M2Q1[$count]}` \$ \&"
            mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                mon+=" \& ";
            else
                mon+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

    done
done

count=0
for (( i=0; i<=3; i++ )) ; do
    m=$((i+1))
    for (( j=$m; j<=3; j++ )) ; do
        n=$((j+1))

        x=${M2C[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            mon+="\$ t_{$n$n$m} \$ \&"
            mon+=" \$ `round 5 ${M2C[$count]}` \$ \&"
            mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                mon+=" \& ";
            else
                mon+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

        x=${M2C[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            mon+="\$ t_{$m$m$n} \$ \&"
            mon+=" \$ `round 5 ${M2C[$count]}` \$ \&"
            mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                mon+=" \& ";
            else
                mon+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

    done
done

count=0
for (( i=0; i<=3; i++ )) ; do
    m=$((i+1))
    for (( j=$m; j<=3; j++ )) ; do
        n=$((j+1))

        x=${M2Q2[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            mon+="\$ u_{$n$n$n$m} \$ \&"
            mon+=" \$ `round 5 ${M2Q2[$count]}` \$ \&"
            mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                mon+=" \& ";
            else
                mon+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

        x=${M2Q2[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            mon+="\$ u_{$m$m$m$n} \$ \&"
            mon+=" \$ `round 5 ${M2Q2[$count]}` \$ \&"
            mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                mon+=" \& ";
            else
                mon+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

    done
done

#3MR
count=0
for (( i=0; i<4; i++ )) ; do 
    l=$((i+1))
    for (( j=$l; j<4; j++ )) ; do
        m=$((j+1))
        for (( k=$m; k<4; k++ )) ; do
            n=$((k+1))

            x=${M3C[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                mon+="\$ t_{$l$m$n} \$ \&"
                mon+=" \$ `round 5 ${M3C[$count]}` \$ \&"
                mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    mon+=" \& ";
                else
                    mon+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))
        done
    done
done

count=0
for (( i=0; i<=3; i++ )) ; do 
    l=$((i+1))
    for (( j=$l; j<=3; j++ )) ; do
        m=$((j+1))
        for (( k=$m; k<=3; k++ )) ; do
            n=$((k+1))

            x=${M3Q[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                mon+="\$ u_{$n$n$m$l} \$ \&"
                mon+=" \$ `round 5 ${M3Q[$count]}` \$ \&"
                mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    mon+=" \& ";
                else
                    mon+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))

            x=${M3Q[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                mon+="\$ u_{$m$m$l$n} \$ \&"
                mon+=" \$ `round 5 ${M3Q[$count]}` \$ \&"
                mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    mon+=" \& ";
                else
                    mon+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))

            x=${M3Q[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                mon+="\$ u_{$l$l$n$m} \$ \&"
                mon+=" \$ `round 5 ${M3Q[$count]}` \$ \&"
                mon+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    mon+=" \& ";
                else
                    mon+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))

        done
    done
done

#1MR

col=0
count=0
d1mr=""
for i in {0..7} ; do
    n=$((i+1))

    x=${D1H[$i]}
    if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
        d1mr+="\$ h_{$n$n} \$ \&"
        d1mr+=" \$ `round 5 ${D1H[$i]}` \$ \&"
        d1mr+=" \$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-2}\\\\ {\\\\rm u}^{-1} \$ "
        if [ $((col%2)) -eq 0 ] ; then 
            d1mr+=" \& ";
        else
            d1mr+=" \\\\\\\\ ";
        fi
        col=$((col+1))
    fi
done

for i in {0..7} ; do
    n=$((i+1))

    x=${D1C[$i]}
    if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
        d1mr+="\$ t_{$n$n$n} \$ \&"
        d1mr+=" \$ `round 5 ${D1C[$i]}` \$ \&"
        d1mr+=" \$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
        if [ $((col%2)) -eq 0 ] ; then 
            d1mr+=" \& ";
        else
            d1mr+=" \\\\\\\\ ";
        fi
        col=$((col+1))
    fi
done

for i in {0..7} ; do
    n=$((i+1))

    x=${D1Q[$i]}
    if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
        d1mr+="\$ u_{$n$n$n$n} \$ \&"
        d1mr+=" \$ `round 5 ${D1Q[$i]}` \$ \&"
        d1mr+=" \$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
        if [ $((col%2)) -eq 0 ] ; then 
            d1mr+=" \& ";
        else
            d1mr+=" \\\\\\\\ ";
        fi
        col=$((col+1))
    fi

done

#2MR

d2mrH=""
count=0
col=0
for i in {0..7} ; do 
    m=$((i+1))
    for (( j=$m; j<=7; j++ )) ; do
        n=$((j+1))

        x=${D2H[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            d2mrH+="\$ h_{$n$m} \$ \&"
            d2mrH+=" \$ `round 5 ${D2H[$count]}` \$ \&"
            #d2mrH+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-2}\\\\ {\\\\rm u}^{-1} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                d2mrH+=" \& ";
            else
                d2mrH+=" \\\\\\\\ ";
            fi
            col=0
        fi
        count=$((count+1))

    done
done

d2mrQ1=""
count=0
col=0
for i in {0..7} ; do 
    m=$((i+1))
    for (( j=$m; j<=7; j++ )) ; do
        n=$((j+1))

        x=${D2Q1[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            d2mrQ1+="\$ u_{$n$n$m$m} \$ \&"
            d2mrQ1+=" \$ `round 5 ${D2Q1[$count]}` \$ \&"
            #d2mrQ1+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                d2mrQ1+=" \& ";
            else
                d2mrQ1+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

    done
done


d2mrC=""
count=0
col=0
for (( i=0; i<=7; i++ )) ; do
    m=$((i+1))
    for (( j=$m; j<=7; j++ )) ; do
        n=$((j+1))

        x=${D2C[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            d2mrC+="\$ t_{$n$n$m} \$ \&"
            d2mrC+=" \$ `round 5 ${D2C[$count]}` \$ \&"
            #d2mrC+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                d2mrC+=" \& ";
            else
                d2mrC+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

        x=${D2C[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            d2mrC+="\$ t_{$m$m$n} \$ \&"
            d2mrC+=" \$ `round 5 ${D2C[$count]}` \$ \&"
            #d2mrC+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                d2mrC+=" \& ";
            else
                d2mrC+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

    done
done

d2mrQ2=""
count=0
col=0
for (( i=0; i<=7; i++ )) ; do
    m=$((i+1))
    for (( j=$m; j<=7; j++ )) ; do
        n=$((j+1))

        x=${D2Q2[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            d2mrQ2+="\$ u_{$n$n$n$m} \$ \&"
            d2mrQ2+=" \$ `round 5 ${D2Q2[$count]}` \$ \&"
            #d2mrQ2+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                d2mrQ2+=" \& ";
            else
                d2mrQ2+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

        x=${D2Q2[$count]}
        if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
            d2mrQ2+="\$ u_{$m$m$m$n} \$ \&"
            d2mrQ2+=" \$ `round 5 ${D2Q2[$count]}` \$ \&"
            #d2mrQ2+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
            if [ $((col%2)) -eq 0 ] ; then 
                d2mrQ2+=" \& ";
            else
                d2mrQ2+=" \\\\\\\\ ";
            fi
            col=$((col+1))
        fi
        count=$((count+1))

    done
done

#3MR
d3mrC=""
count=0
col=0
for (( i=0; i<=7; i++ )) ; do 
    l=$((i+1))
    for (( j=$l; j<=7; j++ )) ; do
        m=$((j+1))
        for (( k=$m; k<=7; k++ )) ; do
            n=$((k+1))

            x=${D3C[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                d3mrC+="\$ t_{$l$m$n} \$ \&"
                d3mrC+=" \$ `round 5 ${D3C[$count]}` \$ \&"
                #d3mrC+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-3}\\\\ {\\\\rm u}^{-3/2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    d3mrC+=" \& ";
                else
                    d3mrC+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))
        done
    done
done

d3mrQ=""
count=0
col=0
for (( i=0; i<=7; i++ )) ; do 
    l=$((i+1))
    for (( j=$l; j<=7; j++ )) ; do
        m=$((j+1))
        for (( k=$m; k<=7; k++ )) ; do
            n=$((k+1))

            x=${D3Q[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                d3mrQ+="\$ u_{$n$n$m$l} \$ \&"
                d3mrQ+=" \$ `round 5 ${D3Q[$count]}` \$ \&"
                #d3mrQ+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    d3mrQ+=" \& ";
                else
                    d3mrQ+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))

            x=${D3Q[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                d3mrQ+="\$ u_{$m$m$l$n} \$ \&"
                d3mrQ+=" \$ `round 5 ${D3Q[$count]}` \$ \&"
                #d3mrQ+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    d3mrQ+=" \& ";
                else
                    d3mrQ+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))

            x=${D3Q[$count]}
            if (( $(echo "${x#-} > $THRESHOLD" |bc -l) )); then
                d3mrQ+="\$ u_{$l$l$n$m} \$ \&"
                d3mrQ+=" \$ `round 5 ${D3Q[$count]}` \$ \&"
                #d3mrQ+="\$ E_{\\\\rm h}\\\\ {\\\\rm \\\\AA}^{-4}\\\\ {\\\\rm u}^{-2} \$ "
                if [ $((col%2)) -eq 0 ] ; then 
                    d3mrQ+=" \& ";
                else
                    d3mrQ+=" \\\\\\\\ ";
                fi
                col=$((col+1))
            fi
            count=$((count+1))

        done
    done
done


sed -i '' "s,\?monomer-QFF\?,${mon}," supp/info.tex
sed -i '' "s,\?dimer-QFF-1MR\?,${d1mr}," supp/info.tex
sed -i '' "s,\?dimer-QFF-2MR-Q1\?,${d2mrQ1}," supp/info.tex
sed -i '' "s,\?dimer-QFF-2MR-C\?,${d2mrC}," supp/info.tex
sed -i '' "s,\?dimer-QFF-2MR-Q2\?,${d2mrQ2}," supp/info.tex
sed -i '' "s,\?dimer-QFF-3MR-C\?,${d3mrC}," supp/info.tex
sed -i '' "s,\?dimer-QFF-3MR-Q\?,${d3mrQ}," supp/info.tex

