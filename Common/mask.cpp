#include "mask.h"

// Specialization of the template, which allows for quicker move semantics in the case
// that T is Image::pixelT
template <>
MaskResult<Image::pixelT>::operator Image() && {
	Image re(_rows, _cols, 255, _data);

	_rows = _cols = min = max = 0;
	_data                     = nullptr;

	return re;
}