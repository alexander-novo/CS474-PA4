#pragma once

#include <array>
#include <functional>

#include "image.h"

template <typename T>
class MaskResult {
public:
	MaskResult(unsigned rows, unsigned cols, T min = 0, T max = 0);
	MaskResult(const MaskResult<T>&);  // Copy constructor
	MaskResult(MaskResult<T>&&);       // Move constructor
	~MaskResult();

	operator Image() const&;
	operator Image() &&;
	T* operator[](unsigned i);
	const T* operator[](unsigned i) const;

	/**
	 * Apply a transformation to each pixel.
	 *
	 * @param transformation The transformation to apply.
	 * @param ivp            Whether or not the transformation attains a maximum on
	 *                       the interval [min, max] exactly at one of the end points
	 *                       and a minimum at the other end point. All monotonic
	 *                       transformations have this property.
	 */
	MaskResult<T>& transform(std::function<T(const T&)> transformation,
	                         bool ivp = true);

	T* const& data       = _data;
	const unsigned& rows = _rows;
	const unsigned& cols = _cols;

	T min, max;

private:
	MaskResult(unsigned rows, unsigned cols, T min, T max, T* data);
	T* _data;
	unsigned _rows, _cols;

	void calcMinMax(T& min, T& max) const;
};

template <typename T, std::size_t N>
class Mask {
	static_assert(N % 2 == 1, "Size of mask must be odd!");

public:
	Mask();
	Mask(const T (&values)[N][N]);

	// Convolution
	MaskResult<T> operator*(const Image& image) const;

	T values[N][N];
	// A sum of all the positive and negative values in the mask, respectively
	// Used for remapping to [0,255] after convolution
	const T& posSum = _posSum;
	const T& negSum = _negSum;

private:
	T _posSum;
	T _negSum;
};

template <typename T>
MaskResult<T>::MaskResult(unsigned rows, unsigned cols, T min, T max, T* data)
    : _rows(rows), _cols(cols), min(min), max(max), _data(data) {}

template <typename T>
MaskResult<T>::MaskResult(unsigned rows, unsigned cols, T min, T max)
    : MaskResult(rows, cols, min, max, new T[rows * cols]) {}

template <typename T>
MaskResult<T>::MaskResult(const MaskResult<T>& other)
    : MaskResult(other._rows, other._cols, other.min, other.max) {
#pragma omp parallel for reduction(min : min) reduction(max : max)
	for (unsigned i = 0; i < _rows * _cols; i++) {
		_data[i] = other._data[i];
		if (data[i] < min)
			min = data[i];
		else if (data[i] > max)
			max = data[i];
	}
}

template <typename T>
MaskResult<T>::MaskResult(MaskResult<T>&& other)
    : MaskResult(other._rows, other._cols, other.min, other.max, other._data) {
	other._rows = other._cols = other.min = other.max = 0;
	delete[] other._data;
	other._data = nullptr;
}

template <typename T>
MaskResult<T>::~MaskResult() {
	if (_data != nullptr) delete[] _data;
}

template <typename T>
MaskResult<T>::operator Image() const& {
	Image re(_rows, _cols, 255);

	T min = MaskResult<T>::min, max = MaskResult<T>::max;
	if (min >= max) calcMinMax(min, max);

#pragma omp parallel for
	for (unsigned i = 0; i < _rows * _cols; i++) {
		re.pixels[i] = (data[i] - min) * re.maxVal / (max - min);
	}

	return re;
}

template <typename T>
MaskResult<T>::operator Image() && {
	// If T isn't Image:pixelT, just copy instead of moving
	// 	Image re(_rows, _cols, 255);

	// 	if (min >= max) calcMinMax(min, max);

	// #pragma omp parallel for
	// 	for (unsigned i = 0; i < _rows * _cols; i++) {
	// 		re.pixels[i] = (data[i] - min) * re.maxVal / (max - min);
	// 	}

	// 	return re;

	return MaskResult<T>::operator Image();
}

// Specialization
template <>
MaskResult<Image::pixelT>::operator Image() &&;

template <typename T>
T* MaskResult<T>::operator[](unsigned i) {
	return _data + i * _cols;
}

template <typename T>
const T* MaskResult<T>::operator[](unsigned i) const {
	return _data + i * _cols;
}

template <typename T>
MaskResult<T>& MaskResult<T>::transform(std::function<T(const T&)> transformation,
                                        bool ivp) {
#pragma omp parallel for
	for (unsigned i = 0; i < _rows * _cols; i++) {
		data[i] = transformation(data[i]);
	}

	if (ivp && min < max) {
		min = transformation(min);
		max = transformation(max);
		if (min > max) { std::swap(min, max); }
	} else {
		calcMinMax(min, max);
	}

	return *this;
}

template <typename T>
void MaskResult<T>::calcMinMax(T& min, T& max) const {
	min = std::min(data[0], data[1]);
	max = std::max(data[0], data[1]);

#pragma omp parallel for reduction(min : min) reduction(max : max)
	for (unsigned i = 2; i < _rows * _cols - 1; i += 2) {
		if (data[i] < data[i + 1]) {
			min = std::min(min, data[i]);
			max = std::max(max, data[i + 1]);
		} else {
			min = std::min(min, data[i + 1]);
			max = std::max(max, data[i]);
		}
	}
}

template <typename T, std::size_t N>
Mask<T, N>::Mask() {}

template <typename T, std::size_t N>
Mask<T, N>::Mask(const T (&values)[N][N]) {
	_posSum = _negSum = 0;
	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = 0; j < N; j++) {
			Mask<T, N>::values[i][j] = values[i][j];
			(values[i][j] < 0 ? _negSum : _posSum) += values[i][j];
		}
	}
}

template <typename T, std::size_t N>
MaskResult<T> Mask<T, N>::operator*(const Image& image) const {
	MaskResult<T> re(image.rows, image.cols, _negSum * image.maxVal,
	                 _posSum * image.maxVal);

	// x and y (and u and v) are interchanged from normal for cache locality purposes
#pragma omp parallel for collapse(2)
	for (unsigned y = 0; y < re.rows; y++) {
		for (unsigned x = 0; x < re.cols; x++) {
			// Keep track of a running sum in a temporary variable instead of in the
			// return image, since we probably have higher precision in T than in
			// pixelT.
			T sum = 0;

			for (unsigned v = 0; v < N; v++) {
				// Clamp the convolution to the image, "extending" the image by
				// its border pixels
				unsigned v_mod = (y + N / 2 < v) ?
				                     0 :
				                     ((y + N / 2 >= v + re.rows) ? (re.rows - 1) :
				                                                   (y - v + N / 2));
				for (unsigned u = 0; u < N; u++) {
					// Clamp u too
					unsigned u_mod = (x + N / 2 < u) ? 0 :
					                                   ((x + N / 2 >= u + re.cols) ?
					                                        (re.cols - 1) :
					                                        (x - u + N / 2));

					sum += values[v][u] * image[v_mod][u_mod];
				}
			}

			re[y][x] = sum;
		}
	}

	return re;
}