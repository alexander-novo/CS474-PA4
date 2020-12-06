# UNR CS 474 Programming Assignment 3

## Building the assignment
A Makefile is provided in the root directory which makes all executables and images required for the homework by default (but not the plots required for Experiment 1). Built executables can be found in their corresponding source folders (`Experiment1/`, `Experiment2/`, etc.) and output images can be found in `out/`.

The only prerequisites for default targets (anything required for the homework) are a working `g++` on the path which supports C++ 14 and OpenMP.

## Running Executables

### Experiment 1
`experiment` - will automatically run and generate needed outputs

### Experiment 2
`rect <image height> <image width> <rectangle height> <rectangle width> <output file>` - Will generate an image of the specified dimensions with a centered rectangle of the specified dimensions. See additional target below.

```
Usage: spectrum <input image> <output image> [options]   (1)
   or: spectrum -h                                       (2)

(1) Take an image file as input, and output the spectrum
    of its fourier transform.
(2) Print this help menu

OPTIONS
  -s  Shift the fourier transform to center its period
      on the center of the image.
  -l  Apply a logarithmic intensity transformation to the
      resulting image, allowing for better visualisation of
      the outer edges of the period.
```

### Experiment 3
```
Usage: reconstruct -p  <input image> <output image>          (1)
   or: reconstruct -m  <input image> <output image>          (2)
   or: reconstruct -ml <input image> <output image>          (3)
   or: reconstruct -h                                        (4)

(1) (phase-only) Take an image file as input, calculate
    Fourier Transform, set magnitude to 0, and then reconstruct
    the image through inverse Fourier Transform.
(2) (magnitude-only) Take an image file as input, calculate
    Fourier Transform, set phase to 0, and then reconstruct
    the image through inverse Fourier Transform.
(3) (magnitude-only with logarithmic inenstiy transformation)
    Same as (2) but a logarithmic intensity transformation
    is applied, allowing to see dark details.
(4) Print this help menu
```

## Building the report
The Makefile does not have a target which generates the report, since different people like to use different engines to generate reports from `.tex` files. I recommend using `pdflatex` from either [TeX Live](https://www.tug.org/texlive/) on Linux or [MiKTeX](https://miktex.org/) on Windows.

The Makefile *does*, however, have a `report` target which will generate all of the neccesary prerequisites of the report (such as any output images and diagrams needed by the report). This target has a few more prerequisites than the default target, however. [Gnuplot](http://www.gnuplot.info/) is used to generate plots for Experiment 1; and [`pnmtopng`](http://netpbm.sourceforge.net/doc/pnmtopng.html) is used to convert `.pgm` images (unreadable by TeX) into `.png` images. These can both be installed on Debian/Ubuntu by using `sudo apt install gnuplot netpbm` and can be found in prebuilt binary form for Windows on their websites (see earlier links).

## Additional targets
In addition to files required by the homework and for the report, the Makefile includes targets for all sorts of output files:

- `out/rect_%.pgm` - Makes an image of a rectangle with specific dimensions e.g. `make out/rect_512_512_32_32.pgm` makes a 512x512 image with a 32x32 rectangle.

- `out/spectrum_rect_%.pgm` - Plots the spectrum of a rectangle with specific dimensions e.g. `make out/spectrum_rect_512_512_32_32.pgm` (see above).

- `out/spectrum_shifted_rect_%.pgm` - Same as above, but applies a transformation before taking the Fourier Transform to shift the spectrum to be centeredo n the image e.g. `make out/spectrum_shifted_rect_512_512_32_32.pgm` (see above).

- `out/spectrum_log_rect_%.pgm` - Same as above, but applied a lograithmic intensity transformation for better contrast e.g. `make out/spectrum_log_rect_512_512_32_32.pgm` (see above).

- `out/%_reconstructed_phase.pgm` - Take an image from the `Images/` directory, apply the Fourier Transform, get rid of the magnitude information, and then reconstruct the iamge with the inverse Fourier Transform e.g. `make out/lenna_reconstructed_phase.pgm`.

- `out/%_reconstructed_mag.pgm` - Same as above, but this time get rid of phase information e.g. `make out/lenna_reconstructed_mag.pgm`.

- `out/%_reconstructed_mag_log.pgm` - Same as above, but apply a logarithmic intensity transformation to the results to improve contrast e.g. `make out/lenna_reconstructed_mag_log.pgm`.