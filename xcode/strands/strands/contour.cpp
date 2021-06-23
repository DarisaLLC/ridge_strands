#include "contour.h"
#include <vector>


contour::contour() { clear(); }

//constructor
contour::contour(int32_t nnum, std::vector<float>& nrow, std::vector<float>& ncol, std::vector<float>& nangle, std::vector<float>& nresponse, contour_class ncont_class) {
	num = nnum;
	row = nrow;
	col = ncol;
	angle = nangle;
	response = nresponse;
	cont_class = ncont_class;
}

void contour::clear() {
	num = 0;
	row.resize(0);
	col.resize(0);
	angle.resize(0);
	response.resize(0);
	width_l.resize(0);
	width_r.resize(0);
	asymmetry.resize(0);
	contrast.resize(0);
	cont_class = contour_class::cont_no_junc;
}
