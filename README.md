# NMR

This project stores my code for the CSC 335 NMR project.

## Compiling

Type `make`, and the program will compile itself. The output will be `nmr.out`.

## Running

The program will read in the program options from standard input using the
format specified in the project description. It is recommended to use file
redirection (like shown below) instead of specifying the options from
the keyboard. Results will be printed to the output file specified in the
program options. Any errors will be printed to standard output.

Example Usage: `nmr.out < nmr.in`.

## Program Options
These are read from standard input (or from a file if file redirection is
used.) They are specified in the following order with one entry on each
line:

1. input file: This is the name of the file that contains the raw spectrum
of points. Each point should be specified in terms of its x and y components
with one point per line.

2. baseline adjustment: Peaks are determined based on points having a y-value
above this real number.

3. tolerance: This real number specifies the tolerance for numerical
algorithms where applicable.

4. filter type: This integer indicates what kind of filtering should be
applied to the raw data. The options are: 0 (no filter), 1 (boxcar filter),
2 (Savitzsky-Golay filter), and 3 (discrete Fourier transform filter).

5. filter size or recovery method: This should be 0 if no filter is selected.
If the boxcar filter is selected, this number should be odd and indicates the
number of points used when averaging. If the Savitzsky-Golay filter is
selected, this number should be one of 5, 11, or 17, and it specifies the
number of points used when averaging. If the discrete Fourier transform filter
is selected, this number indicates the way that the filtered points are
recovered from the Fourier domain and should be one of 0 (inverse method),
1 (direct method via Gaussian elimination), or 2 (direct method via the
Jacobi iterative method).

6. filter passes: This should be 0 if no filter or the discrete Fourier
transform filter is selected. Otherwise, this integer indicates the number
of times the filter is applied to the input.

7. integration method: This integer indicates the way that the area of each
peak is computed, and it should be one of the following: 0 (adaptive
quadrature), 1 (Romberg integration), 2 (composite Newton-Cotes), or
3 (Gaussian quadrature).

8. output file: This is the name of the file that the analysis results should
be printed to.

