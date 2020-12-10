// Experiment3/main.cpp
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "../Common/fft.h"
#include "../Common/image.h"
#include "../Common/mask.h"

// Struct for inputting arguments from command line
struct Arguments {
	char *inputImagePath, *outImagePath;
	Image inputImage;
	std::ofstream outFile;
	float gamma_H = 1.5, gamma_L = .5, c = 1, D_0 = 1.8;
};

bool verifyArguments(int argc, char** argv, Arguments& arg, int& err);
void printHelp();

int main(int argc, char** argv) {
	int err;
	Arguments arg;

	if (!verifyArguments(argc, argv, arg, err)) { return err; }

	MaskResult<float> input = arg.inputImage;

	input.transform([](const float& x) { return log(x); });

	std::complex<float>* four = new std::complex<float>[input.rows * input.cols];
	fft2D(input, four);

	float coef      = -arg.c / arg.D_0 / arg.D_0;
	float gammaDiff = arg.gamma_H - arg.gamma_L;
	// Apply high-pass filtering filter via pointwise multiplication
#pragma omp parallel for collapse(2)
	for (unsigned y = 0; y < input.rows; y++) {
		for (unsigned x = 0; x < input.cols; x++) {
			int x_adj = x - input.cols / 2, y_adj = y - input.rows / 2;
			four[y * input.cols + x] *=
			    gammaDiff * (1 - exp(coef * (x_adj * x_adj + y_adj * y_adj))) +
			    arg.gamma_L;
		}
	}

	// Take inverse transform
	fft2D(four, input.rows, input.cols, 1);
	float min, max;
	min = max = four[0].real();
#pragma omp parallel for collapse(2) reduction(min : min) reduction(max : max)
	for (unsigned y = 0; y < input.rows; y++) {
		for (unsigned x = 0; x < input.cols; x++) {
			float shift = ((x + y) % 2 == 0) ? 1 : -1;
			input[y][x] = four[y * input.cols + x].real() * shift;

			if (input[y][x] < min)
				min = input[y][x];
			else if (input[y][x] > max)
				max = input[y][x];
		}
	}
	input.min = min;
	input.max = max;

	input.transform([](const float& x) { return exp(x); });
	arg.outFile << (Image) input;
}

bool verifyArguments(int argc, char** argv, Arguments& arg, int& err) {
	if (argc < 2 ||
	    (argc < 3 && strcmp(argv[1], "-h") && strcmp(argv[1], "--help"))) {
		std::cout << "Missing operand.\n\n";
		err = 1;
		printHelp();
		return false;
	}

	if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
		printHelp();
		return false;
	}

	// Required arguments
	arg.inputImagePath = argv[1];
	std::ifstream inFile(argv[1]);
	try {
		arg.inputImage = Image::read(inFile);
	} catch (std::exception& e) {
		std::cout << "Image \"" << argv[1] << "\" failed to be read: \"" << e.what()
		          << "\"\n";
		err = 2;
		return false;
	}

	arg.outImagePath = argv[argc - 1];
	arg.outFile.open(argv[argc - 1]);
	if (!arg.outFile) {
		std::cout << "Could not open \"" << argv[argc - 1] << "\"\n";
		err = 2;
		return false;
	}

	// Optional Arguments
	for (unsigned i = 2; i < argc - 1; i++) {
		if (!strcmp(argv[i], "-gl")) {
			if (i + 1 >= argc) {
				std::cout << "Missing gamma_L value\n\n";
				err = 1;
				printHelp();
				return false;
			}

			char* end;
			arg.gamma_L = strtof(argv[i + 1], &end);
			if (end == argv[i + 1]) {
				std::cout << "\"" << argv[i + 1]
				          << "\" could not be interpreted as a number\n";
				err = 2;
				return false;
			} else if (arg.gamma_L < 0 || arg.gamma_L > 1) {
				std::cout << "gamma_L = " << arg.gamma_L
				          << " must be between 0 and 1.\n";
				err = 2;
				return false;
			}

			i++;
		} else if (!strcmp(argv[i], "-gh")) {
			if (i + 1 >= argc) {
				std::cout << "Missing gamma_H value\n\n";
				err = 1;
				printHelp();
				return false;
			}

			char* end;
			arg.gamma_H = strtof(argv[i + 1], &end);
			if (end == argv[i + 1]) {
				std::cout << "\"" << argv[i + 1]
				          << "\" could not be interpreted as a number\n";
				err = 2;
				return false;
			} else if (arg.gamma_H < 1 || arg.gamma_H > 2) {
				std::cout << "gamma_H = " << arg.gamma_H
				          << " must be between 1 and 2.\n";
				err = 2;
				return false;
			}

			i++;
		} else if (!strcmp(argv[i], "-c")) {
			if (i + 1 >= argc) {
				std::cout << "Missing c value\n\n";
				err = 1;
				printHelp();
				return false;
			}

			char* end;
			arg.c = strtof(argv[i + 1], &end);
			if (end == argv[i + 1]) {
				std::cout << "\"" << argv[i + 1]
				          << "\" could not be interpreted as a number\n";
				err = 2;
				return false;
			}

			i++;
		} else if (!strcmp(argv[i], "-d0")) {
			if (i + 1 >= argc) {
				std::cout << "Missing D_0 value\n\n";
				err = 1;
				printHelp();
				return false;
			}

			char* end;
			arg.D_0 = strtof(argv[i + 1], &end);
			if (end == argv[i + 1]) {
				std::cout << "\"" << argv[i + 1]
				          << "\" could not be interpreted as a number\n";
				err = 2;
				return false;
			}

			i++;
		} else {
			std::cout << "Unrecognised argument \"" << argv[i] << "\".\n";
			printHelp();
			err = 1;
			return false;
		}
	}

	return true;
}

void printHelp() {
	// Arguments default values
	Arguments def;
	std::cout
	    << "Usage: homomorphic <input image> [options] <output image>    (1)\n"
	    << "   or: homomorphic -h                                        (4)\n\n"
	    << "(1) Take an image file and apply homomorphic filtering to it.\n"
	    << "    Then write the resulting image to the output file.\n"
	    << "(2) Print this help menu\n\n"
	    << "OPTIONS:\n"
	    << "  -gl <value>    Set the gamma_L value. <value> must be able to be\n"
	    << "                 parsed as a float. Default is " << def.gamma_L << ".\n"
	    << "  -gh <value>    Set the gamma_H value. <value> must be able to be\n"
	    << "                 parsed as a float. Default is " << def.gamma_H << ".\n"
	    << "  -c <value>     Set the c value. <value> must be able to be\n"
	    << "                 parsed as a float. Default is " << def.c << ".\n"
	    << "  -d0 <value>    Set the D_0 value. <value> must be able to be\n"
	    << "                 parsed as a flaot. Default is " << def.D_0 << ".\n";
}