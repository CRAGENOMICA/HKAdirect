cd ./source
gcc *.c -lm -o HKAdirect -lm -Wall -pedantic
cp ./HKAdirect ../bin

./HKAdirect ../examples/checkHKA.txt  > ../examples/checkHKA_results.txt
./HKAdirect ../examples/checkHKA_ChromX.txt > ../examples/checkHKA_ChromX_results.txt
./HKAdirect ../examples/checkHKA_missing.txt > ../examples/checkHKA_missing_results.txt
./HKAdirect ../examples/HKA1000loci.txt > ../examples/HKA1000loci_results.txt

rm ./HKAdirect
cd ..
