



#include "link.h"
#include "utils.hpp"

using namespace std;

/*
 * Extract the line width by using a facet model line detector on an image of
 * the absolute value of the gradient.
 */

void link::compute_line_width(std::vector<float>& dx, std::vector<float>& dy, bool correct_pos, std::vector<contour>& contours){
	auto width = m_width;
	auto height = m_height;
	int i, j, k;
	int r, c, l;
	int x, y, dir;
	float length;
	contour cont;
	int num_points, max_num_points;
	float d, dr, dc, drr, drc, dcc;
	float i1, i2, i3, i4, i5, i6, i7, i8, i9;
	float t1, t2, t3, t4, t5, t6;
	float a, b, t = 0;
	int num = 0;
	float nx, ny;
	float n1, n2;
	float p1, p2;
	float val;
	float px, py;

	max_num_points = 0;
	for (auto i = 0; i < contours.size(); i++) {
		num_points = contours[i].num;
		if (num_points > max_num_points)
			max_num_points = num_points;
	}

	std::vector<float> width_l(max_num_points);
	std::vector<float> width_r(max_num_points);
	std::vector<float> grad_l(max_num_points);
	std::vector<float> grad_r(max_num_points);
	std::vector<float> pos_x(max_num_points);
	std::vector<float> pos_y(max_num_points);
	std::vector<float> correct(max_num_points);
	std::vector<float> contrast(max_num_points);
	std::vector<float> asymm(max_num_points);

	std::vector<float> grad(width * height);
	length = 2.5f * sigma();

#if 0

	max_line = (int)Math.ceil(length * 3);
	line = new Offset[max_line];
	for (int o = 0; o < line.length; o++) {
		line[o] = new Offset();
	}
#endif

	/* Compute the gradient image. */
	for (r = 0; r < height; r++) {
		for (c = 0; c < width; c++) {
			l = LCOR(r, c, width);
			grad[l] = (float)sqrt(dx[l] * dx[l] + dy[l] * dy[l]);
		}
	}

	for (i = 0; i < contours.size(); i++) {
		cont = contours[i];
		num_points = cont.num;

		for (j = 0; j < num_points; j++) {
			px = cont.row[j];
			py = cont.col[j];
			pos_x[j] = px;
			pos_y[j] = py;
			r = (int)floor(px + 0.5);
			c = (int)floor(py + 0.5);
			nx = cos(cont.angle[j]);
			ny = sin(cont.angle[j]);
			/* Compute the search line. */
			std::vector<offset> line;
			bresenham(nx, ny, 0.0, 0.0, length, line);
			auto num_line = line.size();
			width_r[j] = width_l[j] = 0;
			/* Look on both sides of the line. */
			for (dir = -1; dir <= 1; dir += 2) {
				for (k = 0; k < num_line; k++) {
					x = BR(r + dir * line[k].x());
					y = BC(c + dir * line[k].y());
					i1 = grad[LCOR(BR(x - 1), BC(y - 1), width)];
					i2 = grad[LCOR(BR(x - 1), y, width)];
					i3 = grad[LCOR(BR(x - 1), BC(y + 1), width)];
					i4 = grad[LCOR(x, BC(y - 1), width)];
					i5 = grad[LCOR(x, y, width)];
					i6 = grad[LCOR(x, BC(y + 1), width)];
					i7 = grad[LCOR(BR(x + 1), BC(y - 1), width)];
					i8 = grad[LCOR(BR(x + 1), y, width)];
					i9 = grad[LCOR(BR(x + 1), BC(y + 1), width)];
					t1 = i1 + i2 + i3;
					t2 = i4 + i5 + i6;
					t3 = i7 + i8 + i9;
					t4 = i1 + i4 + i7;
					t5 = i2 + i5 + i8;
					t6 = i3 + i6 + i9;
					dr = (t3 - t1) / 6;
					dc = (t6 - t4) / 6;
					drr = (t1 - 2 * t2 + t3) / 6;
					dcc = (t4 - 2 * t5 + t6) / 6;
					drc = (i1 - i3 - i7 + i9) / 4;
					std::vector<std::vector<float>> eigenvalvect = compute_eigenvals(2 * drr, drc, 2 * dcc);
					val = -eigenvalvect[0][0];
					if (val > 0.0) {
						n1 = eigenvalvect[1][0]; // igenvector coordinate y
						n2 = eigenvalvect[1][1]; // eigenvector coordinate x
						a = 2.0 * (drr * n1 * n1 + drc * n1 * n2 + dcc * n2 * n2);
						b = dr * n1 + dc * n2;
						if (a == 0.0) continue;
						t = (-1) * b / a;
						p1 = t * n1;
						p2 = t * n2;
						if (abs(p1) <= 0.5 && abs(p2) <= 0.5) {
							/*
								* Project the maximum point position perpendicularly onto the search line.
								*/
							a = 1;
							b = nx * (px - (r + dir * line[k].x() + p1)) + ny * (py - (c + dir * line[k].y() + p2));
							if (a == 0.0) continue;
							t = (-1) * b / a;
							d = (-i1 + 2 * i2 - i3 + 2 * i4 + 5 * i5 + 2 * i6 - i7 + 2 * i8 - i9) / 9;
							if (dir == 1) {
								grad_r[j] = d + p1 * dr + p2 * dc + p1 * p1 * drr + p1 * p2 * drc
									+ p2 * p2 * dcc;
								width_r[j] = abs(t);
							}
							else {
								grad_l[j] = d + p1 * dr + p2 * dc + p1 * p1 * drr + p1 * p2 * drc
									+ p2 * p2 * dcc;
								width_l[j] = abs(t);
							}
							break;
						}
					}
				}
			}
		}

	//	fix_locations(width_l, width_r, grad_l, grad_r, pos_x, pos_y, correct, contrast, asymm, sigma, mode,
//			correct_pos, cont);
	}
}


/*
 * Fill gaps in the arrays master, slave1, and slave2, i.e., points where
 * master=0, by interpolation (interior points) or extrapolation (end points).
 * The array master will usually be the width of the line, while slave1 and
 * slave2 will be values that depend on master[i] being 0, e.g., the gradient at
 * each line point. The arrays slave1 and slave2 can be NULL.
 */

void link::fill_gaps(std::vector<float>& master, std::vector<float>& slave1, std::vector<float>& slave2, contour& cont){
	int i, j, k, s, e;
	int num_points;
	float m_s, m_e, s1_s, s1_e, s2_s, s2_e, d_r, d_c, arc_len, len;

	num_points = cont.num;
	for (i = 0; i < num_points; i++) {
		if (master[i] == 0) {
			for (j = i + 1; j < num_points; j++) {
				if (master[j] > 0)
					break;
			}
			m_s = 0;
			m_e = 0;
			s1_s = 0;
			s1_e = 0;
			s2_s = 0;
			s2_e = 0;
			if (i > 0 && j < num_points - 1) {
				s = i;
				e = j - 1;
				m_s = master[(s - 1)];
				m_e = master[(e + 1)];
				if (!slave1.empty()) {
					s1_s = slave1[(s - 1)];
					s1_e = slave1[(e + 1)];
				}
				if (!slave2.empty()) {
					s2_s = slave2[(s - 1)];
					s2_e = slave2[(e + 1)];
				}
			}
			else if (i > 0) {
				s = i;
				e = num_points - 2;
				m_s = master[(s - 1)];
				m_e = master[(s - 1)];
				master[(e + 1)] = m_e;
				if (!slave1.empty()) {
					s1_s = slave1[(s - 1)];
					s1_e = slave1[(s - 1)];
					slave1[(e + 1)] = s1_e;
				}
				if (!slave2.empty()) {
					s2_s = slave2[(s - 1)];
					s2_e = slave2[(s - 1)];
					slave2[(e + 1)] = s2_e;
				}
			}
			else if (j < num_points - 1) {
				s = 1;
				e = j - 1;
				m_s = master[(e + 1)];
				m_e = master[(e + 1)];
				master[(s - 1)] = m_s;
				if (!slave1.empty()) {
					s1_s = slave1[(e + 1)];
					s1_e = slave1[(e + 1)];
					slave1[(s - 1)] = s1_s;
				}
				if (!slave2.empty()) {
					s2_s = slave2[(e + 1)];
					s2_e = slave2[(e + 1)];
					slave2[(s - 1)] = s2_s;
				}
			}
			else {
				s = 1;
				e = num_points - 2;
				m_s = master[(s - 1)];
				m_e = master[(e + 1)];
				if (!slave1.empty()) {
					s1_s = slave1[(s - 1)];
					s1_e = slave1[(e + 1)];
				}
				if (!slave2.empty()) {
					s2_s = slave2[(s - 1)];
					s2_e = slave2[(e + 1)];
				}
			}
			arc_len = 0;
			for (k = s; k <= e + 1; k++) {
				d_r = cont.row[k] - cont.row[(k - 1)];
				d_c = cont.col[k] - cont.col[(k - 1)];
				arc_len += sqrt(d_r * d_r + d_c * d_c);
			}
			len = 0;
			for (k = s; k <= e; k++) {
				d_r = cont.row[k] - cont.row[(k - 1)];
				d_c = cont.col[k] - cont.col[(k - 1)];
				len += sqrt(d_r * d_r + d_c * d_c);
				master[k] = (arc_len - len) / arc_len * m_s + len / arc_len * m_e;
				if (!slave1.empty())
					slave1[k] = (arc_len - len) / arc_len * s1_s + len / arc_len * s1_e;
				if (!slave2.empty())
					slave2[k] = (arc_len - len) / arc_len * s2_s + len / arc_len * s2_e;
			}
			i = j;
		}
	}
}


/*
 * Correct the extracted line positions and widths. The algorithm first closes
 * gaps in the extracted data width_l, width_r, grad_l, and grad_r to provide
 * meaningful input over the whole line. Then the correction is calculated.
 * After this, gaps that have been introduced by the width correction are again
 * closed. Finally, the position correction is applied if correct_pos is set.
 * The results are returned in width_l, width_r, and cont.
 */
void link::fix_locations(std::vector<float>& width_l, std::vector<float>& width_r, std::vector<float>& grad_l, std::vector<float>& grad_r, std::vector<float>& pos_x,
	std::vector<float>& pos_y, std::vector<float>& correction, std::vector<float>& contr, std::vector<float>& asymm, bool correct_pos, contour& cont) {

#if 0
	int i;
	int num_points;
	double px, py;
	double nx, ny;
	double w_est, r_est;
	MutableDouble w_real, h_real, corr;
	w_real = new MutableDouble();
	h_real = new MutableDouble();
	corr = new MutableDouble();
	MutableDouble w_strong = new MutableDouble();
	MutableDouble w_weak = new MutableDouble();
	double correct, asymmetry, response, width, contrast;
	bool weak_is_r;
	bool correct_start, correct_end;
	Convol convol = new Convol();
	fill_gaps(width_l, grad_l, null, cont);
	fill_gaps(width_r, grad_r, null, cont);

	num_points = cont.num;

	/* Calculate true line width, asymmetry, and position correction. */
	if (correct_pos) {
		/*
		 * Do not correct the position of a junction point if its width is found by
		 * interpolation, i.e., if the position could be corrected differently for each
		 * junction point, thereby destroying the junction.
		 */
		correct_start = ((cont.cont_class == contour_class::cont_no_junc
			|| cont.cont_class == contour_class::cont_end_junc
			|| cont.cont_class == contour_class::cont_closed)
			&& (width_r[0] > 0 && width_l[0] > 0));
		correct_end = ((cont.cont_class == contour_class::cont_no_junc
			|| cont.cont_class == contour_class::cont_start_junc
			|| cont.cont_class == contour_class::cont_closed)
			&& (width_r[(num_points - 1)] > 0 && width_l[(num_points - 1)] > 0));
		/*
		 * Calculate the true width and assymetry, and its corresponding correction for
		 * each line point.
		 */
		for (i = 0; i < num_points; i++) {
			if (width_r[i] > 0 && width_l[i] > 0) {
				w_est = (width_r[i] + width_l[i]) * LINE_WIDTH_COMPENSATION;
				if (grad_r[i] <= grad_l[i]) {
					r_est = grad_r[i] / grad_l[i];
					weak_is_r = true;
				}
				else {
					r_est = grad_l[i] / grad_r[i];
					weak_is_r = false;
				}
				Correct.line_corrections(sigma, w_est, r_est, w_real, h_real, corr, w_strong, w_weak);
				w_real.setValue(w_real.getValue() / LINE_WIDTH_COMPENSATION);
				corr.setValue(corr.getValue() / LINE_WIDTH_COMPENSATION);
				width_r[i] = w_real.getValue();
				width_l[i] = w_real.getValue();
				if (weak_is_r) {
					asymm[i] = h_real.getValue();
					correction[i] = -corr.getValue();
				}
				else {
					asymm[i] = -h_real.getValue();
					correction[i] = corr.getValue();
				}
			}
		}

		fill_gaps(width_l, correction, asymm, cont);
		for (i = 0; i < num_points; i++)
			width_r[i] = width_l[i];

		/* Adapt the correction for junction points if necessary. */
		if (!correct_start)
			correction[0] = 0;
		if (!correct_end)
			correction[(num_points - 1)] = 0;

		for (i = 0; i < num_points; i++) {
			px = pos_x[i];
			py = pos_y[i];
			nx = Math.cos(cont.angle[i]);
			ny = Math.sin(cont.angle[i]);
			px = px + correction[i] * nx;
			py = py + correction[i] * ny;
			pos_x[i] = px;
			pos_y[i] = py;
		}
	}

	/* Update the position of a line and add the extracted width. */
	cont.width_l = new float[num_points];
	cont.width_r = new float[num_points];
	for (i = 0; i < num_points; i++) {
		cont.width_l[i] = (float)width_l[i];
		cont.width_r[i] = (float)width_r[i];
		cont.row[i] = (float)pos_x[i];
		cont.col[i] = (float)pos_y[i];
	}

	/* Now calculate the true contrast. */
	if (correct_pos) {
		cont.asymmetry = new float[num_points];
		cont.intensity = new float[num_points];
		for (i = 0; i < num_points; i++) {
			response = cont.response[i];
			asymmetry = Math.abs(asymm[i]);
			correct = Math.abs(correction[i]);
			width = cont.width_l[i];
			if (width < MIN_LINE_WIDTH)
				contrast = 0;
			else
				contrast = (response / Math.abs(convol.phi2(correct + width, sigma)
					+ (asymmetry - 1) * convol.phi2(correct - width, sigma)));

			if (contrast > MAX_CONTRAST)
				contrast = 0;
			contr[i] = contrast;
		}
		fill_gaps(contr, null, null, cont);
		for (i = 0; i < num_points; i++) {
			cont.asymmetry[i] = (float)asymm[i];
			if (mode == MODE_LIGHT)
				cont.intensity[i] = (float)contr[i];
			else
				cont.intensity[i] = (float)-contr[i];
		}
	}
#endif

}

