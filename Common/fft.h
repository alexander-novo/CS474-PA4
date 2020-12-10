#include <complex>

#include "image.h"
#include "mask.h"

/**
 * Fast Fourier Transform
 *
 * @param[in,out] data  A pointer to a list of std::complex<float> which stores
 *                      original function values and ends up storing fourier transform
 *                      values.
 * @param[in]     n     Length of data. Must be a power of 2. Pad with zeroes for
 *                      calculation on arrays of length not a power of 2.
 * @param[in]     isign Either -1 (for forward transform) or 1 (for inverse
 *                      transform).
 * @param[in]     byrow The row size for a column-wise transform.
 */
void fft(std::complex<float> data[], unsigned n, int isign, unsigned byrow = 1);

/**
 * Fast Fourier Transform - parallel version. Much faster at 1D DFT than normal fft,
 * but use only for computing 1D DFT, as higher dimensions can be parallelized much
 * better using separability property. Only uses 4 threads instead of highest number
 * available due to high synchronization needed among threads (and synchronization
 * becomes harder and harder to get with more threads).
 *
 * @param[in,out] data  A pointer to a list of std::complex<float> which stores
 *                      original function values and ends up storing fourier transform
 *                      values.
 * @param[in]     n     Length of data. Must be a power of 2. Pad with zeroes for
 *                      calculation on arrays of length not a power of 2.
 * @param[in]     isign Either -1 (for forward transform) or 1 (for inverse
 *                      transform).
 * @param[in]     byrow The row size for a column-wise transform.
 */
void parallelFFT(std::complex<float> data[], unsigned n, int isign,
                 unsigned byrow = 1);

/**
 * 2 Dimensional Fast Fourier Transform. Uses fft() and separability to compute.
 *
 * @param[in,out] data  A pointer to a list of std::complex<float> which stores
 *                      original function values and ends up storing fourier
 * transform values. Stored in row-major order.
 * @param[in]     rows  Number of rows of data. Must be a power of 2
 * @param[in]     cols  Number of columns of data. Must be a power of 2.
 * @param[in]     isign Either -1 (for forward transform) or 1 (for inverse
 *                      transform).
 */
void fft2D(std::complex<float> data[], unsigned rows, unsigned cols, int isign);

/**
 * 2 Dimensional Fast Fourier Transform. Signature for class purposes.
 *
 * @param[in]     N        Number of columns of data. Must be a power of 2.
 * @param[in]     M        Number of rows of data. Must be a power of 2.
 * @param[in,out] real_Fuv The list of real values of the samples of the function, as
 *                         well as real values of fourier transform.
 * @param[in,out] imag_Fuv The list of imaginary values of the samples of the
 *                         function, as well as imaginary values of fourier transform.
 * @param[in]     isign    Either -1 (for forward transform) or 1 (for inverse
 *                         transform).
 */
void fft2D(unsigned N, unsigned M, float real_Fuv[], float imag_Fuv[], int isign);

/**
 * 2 Dimensional Fast Fourier Transform of an image. Only computes forward transform.
 * To take inverse transform, just take inverse normally and calculate image from
 * results.
 *
 * @param[in]  im         An image to take the fourier transform of.
 * @param[out] transform  Where the fourier transform is stored. Must have a capacity
 *                        of at least im.rows * im.cols
 */
void fft2D(const Image& im, std::complex<float> transform[], bool shift = true);

/**
 * 2 Dimensional Fast Fourier Transform of a MaskResult. Only computes forward
 * transform. To take inverse transform, just take inverse normally and calculate
 * image from results.
 *
 * @param[in]  im         A mask result to take the fourier transform of.
 * @param[out] transform  Where the fourier transform is stored. Must have a capacity
 *                        of at least im.rows * im.cols
 */
void fft2D(const MaskResult<float>& im, std::complex<float> transform[],
           bool shift = true);

/**
 * Plot the spectrum of a 2D Fourier Transform
 *
 *
 */
MaskResult<float> plotSpectrum(std::complex<float> transform[], unsigned rows,
                               unsigned cols);