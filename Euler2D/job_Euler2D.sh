mkdir DATA
mkdir animations

g++ -O3 -Wall Euler2D.cpp -o cylindrical_explosion
./cylindrical_explosion

rm cylindrical_explosion

mv *.dat DATA

echo "Initialization of animations"
python3 plot2d.py 
mv *.gif animations
echo "Animations done!