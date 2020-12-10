// Experiment2/main.cpp
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../Common/fft.h"
#include "../Common/image.h"

const float sobel_mask[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};

void normalize_image(Image& image, float* pixels) {
	float max = -1000000.0;
	float min = 10000000.0;

	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			float value = pixels[i * image.rows + j];

			max = std::max(value, max);
			min = std::min(value, min);
		}
	}

	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			image[i][j] =
			    (int) (255.0 * (((double) pixels[i * image.rows + j] - min) /
			                    (double) (max - min)));
		}
	}
}

void spatially_filter(Image& image, int mask_size) {
	float pixels[image.cols * image.rows];

	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			float output_pixel_value = 0;

			for (int k = -1; k < 2; k++) {
				for (int l = -1; l < 2; l++) {
					if (i + k < 0 || i + k >= image.rows || j + l < 0 ||
					    j + l >= image.cols) {
						output_pixel_value += 0;
					}
					// calculate output value
					else
						output_pixel_value +=
						    image[i + k][j + l] *
						    sobel_mask[k + mask_size / 2][l + mask_size / 2];

					// std::cout << k + mask_size/2 << std::endl;
				}
			}

			pixels[i * image.cols + j] = output_pixel_value;
		}
	}

	normalize_image(image, pixels);
}

void freq_filter(Image& image, Image& maskSpectrum,
                 std::complex<float>* mask_transform) {
	float pixels[image.cols * image.rows];

	// take fft of image with shift
	std::complex<float>* transform = new std::complex<float>[image.rows * image.cols];

	fft2D(image, transform, true);

	// intialize h(x,y)
	float mask[image.rows][image.cols];
	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) { mask[i][j] = 0; }
	}

	// place sobel mask in center of h(x,y), offset by one for zero padding
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int x = image.rows / 2 - 1 + i;
			int y = image.cols / 2 - 1 + j;

			mask[x][y] =
			    sobel_mask[i][j] * pow((-1), x + y);  // center spectrum of mask
		}
	}

	// take fft of h(x,y)
	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			mask_transform[i * image.cols + j] = std::complex<float>(mask[i][j], 0);
		}
	}

	fft2D(mask_transform, image.cols, image.rows, -1);

	// Element-wise complex multiplication
	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			int index = i * image.cols + j;

			// set real part to zero, undo cenetering of sobel mask
			mask_transform[index] =
			    std::complex<float>(0, mask_transform[index].imag()) *
			    std::complex<float>(pow(-1, i + j));

			// multiply F(x,y) and H(x,y)
			transform[index] *= mask_transform[index];
		}
	}

	MaskResult<float> spectrum = plotSpectrum(transform, image.rows, image.cols);

	maskSpectrum = spectrum;

	// inverse fft
	fft2D(transform, image.rows, image.cols, 1);

	// update image
	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			pixels[i * image.cols + j] = transform[i * image.cols + j].real() *
			                             pow(-1, i + j);  // undo spectrum centering
		}
	}

	normalize_image(image, pixels);
}

int main(int argc, char** argv) {
	std::cout << "Experiment 2: Frequency Filtering" << std::endl;

	// Load images
	std::ifstream inFile(argv[1]);
	Image image             = Image::read(inFile);
	Image freq_filter_image = Image(image);
	Image sobelSpectrum;
	std::complex<float>* mask_transform =
	    new std::complex<float>[image.rows * image.cols];

	spatially_filter(image, 3);

	freq_filter(freq_filter_image, sobelSpectrum, mask_transform);

	// Save output images
	std::ofstream outFile;

	outFile.open(std::string(argv[2]).substr(0, std::string(argv[2]).find('.')) +
	             "_spatial.pgm");
	outFile << image;
	outFile.close();

	outFile.open(std::string(argv[2]).substr(0, std::string(argv[2]).find('.')) +
	             "_frequency.pgm");
	outFile << freq_filter_image;
	outFile.close();

	outFile.open("out/sobel_spectrum.pgm");
	outFile << sobelSpectrum;
	outFile.close();

	outFile.open("out/sobel_imag.dat");
	outFile << "#  x    y        iz\n";
	for (unsigned x = 0; x < image.rows; x++) {
		for (unsigned y = 0; y < image.cols; y++) {
			outFile << std::setw(4) << x << ' ' << std::setw(4) << y << ' '
			        << std::setw(12) << std::setprecision(9) << std::fixed
			        << mask_transform[y * image.rows + x].imag() << '\n';
		}
		outFile << '\n';
	}
	delete[] mask_transform;
}