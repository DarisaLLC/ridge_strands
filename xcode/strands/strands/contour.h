#pragma once

#include<vector>
#include <cmath>


/** This type holds one extracted line.  The field num contains the number of
points in the line.  The coordinates of the line points are given in the
arrays row and col.  The array angle contains the direction of the normal
to each line point, as measured from the row-axis.  Some people like to
call the col-axis the x-axis and the row-axis the y-axis, and measure the
angle from the x-axis.  To convert the angle into this convention, subtract
PI/2 from the angle and normalize it to be in the interval [0,2*PI).  The
array response contains the response of the operator, i.e., the second
directional derivative in the direction of angle, at each line point.  The
arrays width_l and width_r contain the width information for each line point
if the algorithm was requested to extract it; otherwise they are NULL.  If
the line position and width correction was applied the contents of width_l
and width_r will be identical.  The arrays asymmetry and contrast contain
the true asymmetry and contrast of each line point if the algorithm was
instructed to apply the width and position correction.  Otherwise, they are
set to NULL.  If the asymmetry, i.e., the weaker gradient, is on the right
side of the line, the asymmetry is set to a positive value, while if it is
on the left side it is set to a negative value. */

/** This type determines the class of a line, i.e., whether its end points are
junction points or whether the line forms a closed loop. */
enum class contour_class
{
	/** no end point is a junction */
	cont_no_junc,
	/** only the start point of the line is a junction */
	cont_start_junc,
	/** only the end point of the line is a junction */
	cont_end_junc,
	/** both end points of the line are junctions */
	cont_both_junc,
	/** the contour is closed */
	cont_closed
};

class contour
{
	/** number of points */
public:

	//default constructor
	contour();

	//constructor
	contour(int32_t nnum, std::vector<float>& nrow, std::vector<float>& ncol, std::vector<float>& nangle, std::vector<float>& nresponse, contour_class ncont_class);

	void clear();

	float compute_length() const {
		std::vector<float>::const_iterator row_b = row.begin() + 1;
		std::vector<float>::const_iterator col_b = col.begin() + 1;
		float length = 0.0f;
		for (; row_b < row.end() && col_b < col.end(); row_b++, col_b++) {
			float dr = *(row_b)-*(row_b - 1);
		    float dc = *(col_b)-*(col_b - 1);
			length += std::sqrt(dr * dr + dc * dc);
		}
		return length;
	}

	int32_t num = 0;
	std::vector<float> row;
	/** column coordinates of the line points (X coordinate in ImageJ)  */
	std::vector<float> col;
	/** angle of normal (measured from the row (Y) axis) */
	std::vector<float> angle;
	/** response of line point (second derivative) */
	std::vector<float> response;
	/** width to the left of the line */
	std::vector<float> width_l;
	/** width to the right of the line */
	std::vector<float> width_r;
	/** asymmetry of the line point */
	std::vector<float> asymmetry;
	/** contrast of the line point */
	std::vector<float> contrast;
	/** contour class (e.g., closed, no_junc) */
	contour_class cont_class;

};


class strand_feature : private contour {

	strand_feature();
	strand_feature(const contour&);






};