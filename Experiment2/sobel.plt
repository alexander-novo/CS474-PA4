set terminal pdf enhanced size 8in, 4.8in
set output 'out/sobel.pdf'

# set hidden3d
set dgrid3d 32,32,
set pm3d at bs
set xlabel 'u'
set ylabel 'v'
set zlabel 'i'
set xrange [0:256]
set yrange [0:256]
set zrange [-.0101:.01]
set cbrange [-.0101:.01]
set xtics 0,64,256
set ytics 0,64,256
set ztics -.01,.01,.01
set cbtics -.01,.01,.01

splot 'out/sobel_imag.dat' u 1:2:3 w pm3d notitle