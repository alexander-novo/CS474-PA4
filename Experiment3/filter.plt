if (!exists("outfile")) outfile='plot.pdf'
if (!exists("gh")) gh=1.5
if (!exists("gl")) gl=.5
if (!exists("d0")) d0=1.8
if (!exists("c")) c=1

set terminal pdf enhanced size 8in, 4.8in
set output outfile
set hidden3d 
set isosamples 100
set pm3d at bs
set ztics ('{/Symbol g}_L' gl+.01, '{/Symbol g}_H' gh)
set cbtics ('{/Symbol g}_L' gl+.01, '{/Symbol g}_H' gh)
set contour both
set cntrparam levels discrete (gh - gl)*(1 - exp(-c))+gl
set clabel 'D_0'
set xrange [-5:5]
set yrange [-5:5]
set zrange [gl-.01:gh]

splot (gh-gl)*(1-exp(-c*(x**2+y**2)/d0**2))+gl title '({/Symbol g}_H - {/Symbol g}_L)(1 - exp(-c(x^2 + y^2) / D_0^2)) + {/Symbol g}_L'