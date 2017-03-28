rm testresults data_*
for l in ??_*; do 
    cd $l
    bash perform.sh
    cd ..
done
echo "Done"
cat testresults