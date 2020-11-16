for file in fig*.gnu
do
    gnuplot -c $file
done

for file in Is-2D_*.eps
do
    epstopdf --autorotate=All $file
done
rm -f *.eps
