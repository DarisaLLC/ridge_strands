#pragma once
#pragma once

#include <vector>
#include <cmath>
#include "strands.h"

class convol
{

	/** Derivative in row direction */
public:

	convol(int32_t cwidth, int32_t cm_height);

	void get_all_derivatives(std::vector<double>& image, double sigma, std::vector<std::vector<float>>& out);

	bool debug() const { return m_debug; }
	void debug(bool state) const { m_debug = state; }


private:

	mutable bool m_debug;
	static constexpr int DERIV_R = 1;
	/** Derivative in column direction */
	static constexpr int DERIV_C = 2;
	/** Second derivative in row direction */
	static constexpr int DERIV_RR = 3;
	/** Second derivative in row and column direction */
	static constexpr int DERIV_RC = 4;
	/** Second derivative in column direction */
	static constexpr int DERIV_CC = 5;

	static double SQRT_2_PI_INV;

	/* Compute the integral of the Gaussian, i.e., the normal distribution. */

	static double SQRTPI;

	static double UPPERLIMIT;

	static double SQRT2;
	static double P10;
	static double P11;
	static double P12;
	static double P13;
	static double Q10;
	static double Q11;
	static double Q12;
	static double Q13;

	static double P20;
	static double P21;
	static double P22;
	static double P23;
	static double P24;
	static double P25;
	static double P26;
	static double P27;
	static double Q20;
	static double Q21;
	static double Q22;
	static double Q23;
	static double Q24;
	static double Q25;
	static double Q26;
	static double Q27;

	static double P30;
	static double P31;
	static double P32;
	static double P33;
	static double P34;
	static double Q30;
	static double Q31;
	static double Q32;
	static double Q33;
	static double Q34;

	/** Size for Gaussian mask */
	static double MAX_SIZE_MASK_0;
	/** Size for 1st derivative mask */
	static double MAX_SIZE_MASK_1;
	/** Size for 2nd derivative mask */
	static double MAX_SIZE_MASK_2;

	static int MASK_SIZE(double MAX, double sigma);



	int32_t width = 0;
	int32_t height = 0;



	static double normal(double x);



	/** Integral of the Gaussian function */
	static double phi0(double x, double sigma);

	/** The Gaussian function */
	static double phi1(double x, double sigma);

	/** First derivative of the Gaussian function */
	static double phi2(double x, double sigma);

	/** Functions to compute the one-dimensional convolution masks of the 0th, 1st,
		and 2nd derivative of the Gaussian kernel for a certain smoothing level
		given by sigma.  The mask is allocated by the function and given as the
		return value.  The caller must ensure that this memory is freed.  The
		output is intended to be used as an array with range [-num:num].  Therefore,
		the caller should add num to the return value.  Examples for the calling
		sequence can be found in convolve_gauss.  Examples for the usage of the
		masks are given in convolve_rows_gauss and convolve_cols_gauss. */

		/* Mask sizes in convol.c */


		/** Gaussian smoothing mask */
	static void compute_gauss_mask_0(double sigma, std::vector<double>& out);

	/** First derivative of Gaussian smoothing mask */
	static void compute_gauss_mask_1(double sigma, std::vector<double>& out);

	/** Second derivative of Gaussian smoothing mask */
	static void compute_gauss_mask_2(double sigma, std::vector<double>& out);


	/** Convolve an image with the derivatives of a Gaussian smoothing kernel.
		Since all of the masks are separable, this is done in two steps in the
		function convolve_gauss.  Firstly, the rows of the image are convolved by
		an appropriate one-dimensional mask in convolve_rows_gauss, yielding an
		intermediate float-image h.  Then the columns of this image are convolved
		by another appropriate mask in convolve_cols_gauss to yield the final
		result k.  At the border of the image the gray values are mirrored. */


		/** Convolve the rows of an image with the derivatives of a Gaussian. **/
	void convolve_rows_gauss(std::vector<double>& image, std::vector<double>& mask, std::vector<float>& out);


	/** Convolve the columns of an image with the derivatives of a Gaussian. */
	void convolve_cols_gauss(std::vector<float>& h, std::vector<double>& mask, std::vector<float>& out);


	/** Convolve an image with a derivative of the Gaussian. */
	void convolve_gauss(std::vector<double>& image, double sigma, int deriv_type, std::vector<float>& out);
	void convolve_gauss_opencv(std::vector<double>& image, double sigma, int deriv_type, std::vector<float>& out);





	/** Mirror the row coordinate at the borders of the image; m_height must be a
		defined variable in the calling function containing the image m_height. */
	int32_t BR(int32_t row);

	/** Mirror the column coordinate at the borders of the image; width must be a
		defined variable in the calling function containing the image width. */
	int32_t BC(int32_t col);

};

