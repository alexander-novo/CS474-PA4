# UNR CS 474 Programming Assignment 4

## Building the assignment
A Makefile is provided in the root directory which makes all executables and images required for the homework by default (but not the plots required for Experiment 1). Built executables can be found in their corresponding source folders (`Experiment1/`, `Experiment2/`, etc.) and output images can be found in `out/`.

The only prerequisites for default targets (anything required for the homework) are a working `g++` on the path which supports C++ 14 and OpenMP.

## Running Executables

### Experiment 1

### Experiment 2

### Experiment 3


## Building the report
The Makefile does not have a target which generates the report, since different people like to use different engines to generate reports from `.tex` files. I recommend using `pdflatex` from either [TeX Live](https://www.tug.org/texlive/) on Linux or [MiKTeX](https://miktex.org/) on Windows.

The Makefile *does*, however, have a `report` target which will generate all of the neccesary prerequisites of the report (such as any output images and diagrams needed by the report). This target has a few more prerequisites than the default target, however. [Gnuplot](http://www.gnuplot.info/) is used to generate plots for Experiment 1; and [`pnmtopng`](http://netpbm.sourceforge.net/doc/pnmtopng.html) is used to convert `.pgm` images (unreadable by TeX) into `.png` images. These can both be installed on Debian/Ubuntu by using `sudo apt install gnuplot netpbm` and can be found in prebuilt binary form for Windows on their websites (see earlier links).

## Additional targets
In addition to files required by the homework and for the report, the Makefile includes targets for all sorts of output files: