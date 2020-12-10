// Experiment1/main.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "../Common/image.h"
#include "../Common/fft.h"

// 7x7 Gaussian mask
static const int mask_7x7[7][7] = {

		{1, 1, 2, 2, 2, 1, 1},
		{1, 2, 2, 4, 2, 2, 1},
		{2, 2, 4, 8, 4, 2, 2},
		{2, 4, 8, 16, 8, 4, 2},
		{2, 2, 4, 8, 4, 2, 2},
		{1, 2, 2, 4, 2, 2, 1},
		{1, 1, 2, 2, 2, 1, 1}
							
						};

//15x15 Gaussian mask
static const int mask_15x15[15][15] = {

	{2, 2,  3,  4,  5,  5,  6,  6,  6,  5,  5,  4,  3, 2, 2},
	{2, 3,  4,  5,  7,  7,  8,  8,  8,  7,  7,  5,  4, 3, 2},
	{3, 4,  6,  7,  9, 10, 10, 11, 10, 10,  9,  7,  6, 4, 3},
	{4, 5,  7,  9, 10, 12, 13, 13, 13, 12, 10,  9,  7, 5, 4},
	{5, 7,  9, 11, 13, 14, 15, 16, 15, 14, 13, 11,  9, 7, 5},
	{5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5},
	{6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6},
	{6, 8, 11, 13, 16, 18, 19, 20, 19, 18, 16, 13, 11, 8, 6},
	{6, 8, 10, 13, 15, 17, 19, 19, 19, 17, 15, 13, 10, 8, 6},
	{5, 7, 10, 12, 14, 16, 17, 18, 17, 16, 14, 12, 10, 7, 5},
	{5, 7,  9, 11, 13, 14, 15, 16, 15, 14, 13, 11,  9, 7, 5},
	{4, 5,  7,  9, 10, 12, 13, 13, 13, 12, 10,  9,  7, 5, 4},
	{3, 4,  6,  7,  9, 10, 10, 11, 10, 10,  9,  7,  6, 4, 3},
	{2, 3,  4,  5,  7,  7,  8,  8,  8,  7,  7,  5,  4, 3, 2},
	{2, 2,  3,  4,  5,  5,  6,  6,  6,  5,  5,  4,  3, 2, 2},
};

void normalize_image(Image& image, float* pixels) {

	float max = -1000000.0;
	float min = 10000000.0;

	for (int i = 0; i < image.rows; i++){
		for (int j = 0; j < image.cols; j++){
		
			float value = pixels[i * image.rows + j];

			max = std::max(value, max);
			min = std::min(value, min);
		}
	}

	for (int i = 0; i < image.rows; i++){
		for (int j = 0; j < image.cols; j++){
		
			image[i][j] = (int) (255.0 * (((double)pixels[i * image.rows + j] - min) / (double)(max - min)));
		}
	}
}

/*
  Smoothes an image based on the Gaussian mask
  @Param: image - the input image that will be smoothed
  @Param: mask_size - the width and height of the mask
  @Return: void
*/
void smooth_image_gaussian(Image& image, int mask_size){

	int normalizion_factor = 0;

	Image originalImage = Image(image);

	// calculate normalizion factor
	for(int i = 0; i < mask_size; i++){
		for (int j = 0; j < mask_size; j++){

			if(mask_size == 7)
				normalizion_factor += mask_7x7[i][j];
			else
				normalizion_factor += mask_15x15[i][j];
		}
	}

	// iterate through image pixels
    for(int i = 0; i < image.cols; i++){
   		for(int j = 0; j < image.rows; j++) {

   			int output_pixel_value = 0;

   			// iterate through mask
   			for (int k = -mask_size/2; k < mask_size/2 + 1; k++)
   			{
   				for (int l = -mask_size/2; l < mask_size/2 + 1; l++)
   				{
   					//bounds checking and padding
   					if(i + k < 0 || i + k >= image.cols || j + l < 0 || j + l >= image.rows){
   						output_pixel_value += 0;
   					}
   					//calculate output value
   					else if(mask_size == 7)
   						output_pixel_value += originalImage[i + k][j + l] * mask_7x7[k + mask_size/2][l + mask_size/2];
   					else
   						output_pixel_value += originalImage[i + k][j + l] * mask_15x15[k + mask_size/2][l + mask_size/2];
   				}
   			}

   			// update image pixel
   			image[i][j] = (int) (output_pixel_value / normalizion_factor);
   		}
   	}
}

/*
  Filters an image in order to remove periodic noise
  @Param: image - the noisy input image that will be filtered
  @Param: return_noise - whether to return just the noise from the image
  @Return: void
*/
void filter_noise(Image& image, bool return_noise) {

	// take FFT of image
	std::complex<float>* transform =
	    new std::complex<float>[image.rows * image.cols];

	fft2D(image, transform, true);

	float pixels[image.rows * image.cols];

	// filter image frequencies
	for (unsigned i = 0; i < image.rows; i++) {
		for (unsigned j = 0; j < image.cols; j++) {

			if(return_noise == false){

				// remove noisy frequencies
				if(abs(i - 256) == 16 && abs(j - 256) == 32)
					transform[i*image.cols + j] = std::complex<float>(0, 0);
			}
			else {

				// remove all frequencies except noisy frequencies
				if( !(abs(i - 256) == 16 && abs(j - 256) == 32))
					transform[i*image.cols + j] = std::complex<float>(0, 0);
			}

			double mag = std::abs(transform[i*image.cols + j]);
			pixels[i*image.rows + j] = 20 * log(1 + mag);
		}
	}

	// inverse fft
	fft2D(transform, image.rows, image.cols, 1);

	// reconstruct image
	for (unsigned i = 0; i < image.rows; i++) {
		for (unsigned j = 0; j < image.cols; j++) {

			pixels[i*image.rows + j] = transform[i*image.rows + j].real() * pow(-1, i + j);
		}
	}

	normalize_image(image, pixels);
}

int main(int argc, char** argv) {
	std::cout << "Experiment 1: Noise Removal" << std::endl;

	// Read and copy input image
	std::ifstream inFile(argv[1]);
 	Image image = Image::read(inFile);
 	Image spectrum = Image(image);
 	Image smoothed_7 = Image(image);
 	Image smoothed_15 = Image(image);
 	Image noise = Image(image);

 	// get spectrum of image for visualization
 	std::complex<float>* transform =
	    new std::complex<float>[image.rows * image.cols];

	fft2D(spectrum, transform, true);

	float pixels[image.rows * image.cols];
	for (unsigned i = 0; i < spectrum.rows; i++) {
		for (unsigned j = 0; j < spectrum.cols; j++) {

			double mag = std::abs(transform[i*image.cols + j]);
			pixels[i*image.rows + j] = 20 * log(1 + mag);
		}
	}
	normalize_image(spectrum, pixels);

	// Filter image
 	filter_noise(image, false);

 	// Return noise from image
 	filter_noise(noise, true);

 	// Smooth orignal image for comparison
 	smooth_image_gaussian(smoothed_7, 7);
	smooth_image_gaussian(smoothed_15, 15);

	// Save output images
	std::ofstream outFile;
	outFile.open(argv[2]);
	outFile << image;
	outFile.close();

	outFile.open(std::string(argv[2]).substr(0, std::string(argv[2]).find('.')) + "_smoothed_7.pgm");
	outFile << smoothed_7;
	outFile.close();

	outFile.open(std::string(argv[2]).substr(0, std::string(argv[2]).find('.')) + "_smoothed_15.pgm");
	outFile << smoothed_15;
	outFile.close();

	outFile.open(std::string(argv[2]).substr(0, std::string(argv[2]).find('.')) + "_noise.pgm");
	outFile << noise;
	outFile.close();

	outFile.open(std::string(argv[2]).substr(0, std::string(argv[2]).find('.')) + "_spectrum.pgm");
	outFile << spectrum;
	outFile.close();

	return 0;
}