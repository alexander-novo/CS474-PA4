set terminal pdf enhanced size 8in, 4.8in
set output 'out/sobel.pdf'

# set hidden3d
set dgrid3d 128,128,
set pm3d at bs
set xlabel 'u'
set ylabel 'v'
set zlabel '𝑖z'
set xrange [1:254]
set yrange [1:254]
set zrange [-.0101:.01]
set cbrange [-.0101:.01]
set xtics 0,32,256
set ytics 0,32,256
set ztics -.01,.01,.01
set cbtics -.01,.01,.01

splot 'out/sobel_imag.dat' u 1:2:3 w pm3d notitle