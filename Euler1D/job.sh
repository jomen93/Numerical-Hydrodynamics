mkdir DATA
g++ -O3 -Wall Euler1.cpp -o Euler1
./Euler1
rm Euler1
mv *.dat DATA
python3 plot.py

