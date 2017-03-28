if [ -d tmp ]; then
    rm -r tmp
fi

NAME=$(pwd | sed s:/:\ :g | awk '{print $(NF)}')

mkdir tmp
cp * tmp
cd tmp
ln -s ../../abinit
/usr/bin/time --output=../../timing_$NAME ./abinit < t1x.files > stdout 
grep etotal oep.00.dat  | tail -n 1 | awk '{print $10}' > result
grep etotal oep.00.dat  | awk -v result=$(cat ../result | awk '{print $1}') '{print $10-result}' > ../../data_$NAME
cat stdout > ../../stdout_$NAME
echo -n $NAME": " >> ../../testresults
paste result ../result | awk -f ../../gradetest.awk >> ../../testresults 
cd ..
#rm -r tmp
