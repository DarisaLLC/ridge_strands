#include "derivatives.h"
#include "utils.hpp"
#include <algorithm>
#include <numeric>

using namespace std;



double convol::SQRT_2_PI_INV = 0.398942280401432677939946059935;
double convol::SQRTPI = 1.772453850905516027;
double convol::UPPERLIMIT = 20.0;
double convol::SQRT2 = 1.41421356237309504880;
double convol::P10 = 242.66795523053175;
double convol::P11 = 21.979261618294152;
double convol::P12 = 6.9963834886191355;
double convol::P13 = -0.035609843701815385;
double convol::Q10 = 215.05887586986120;
double convol::Q11 = 91.164905404514901;
double convol::Q12 = 15.082797630407787;
double convol::Q13 = 1.0;
double convol::P20 = 300.4592610201616005;
double convol::P21 = 451.9189537118729422;
double convol::P22 = 339.3208167343436870;
double convol::P23 = 152.9892850469404039;
double convol::P24 = 43.16222722205673530;
double convol::P25 = 7.211758250883093659;
double convol::P26 = 0.5641955174789739711;
double convol::P27 = -0.0000001368648573827167067;
double convol::Q20 = 300.4592609569832933;
double convol::Q21 = 790.9509253278980272;
double convol::Q22 = 931.3540948506096211;
double convol::Q23 = 638.9802644656311665;
double convol::Q24 = 277.5854447439876434;
double convol::Q25 = 77.00015293522947295;
double convol::Q26 = 12.78272731962942351;
double convol::Q27 = 1.0;
double convol::P30 = -0.00299610707703542174;
double convol::P31 = -0.0494730910623250734;
double convol::P32 = -0.226956593539686930;
double convol::P33 = -0.278661308609647788;
double convol::P34 = -0.0223192459734184686;
double convol::Q30 = 0.0106209230528467918;
double convol::Q31 = 0.191308926107829841;
double convol::Q32 = 1.05167510706793207;
double convol::Q33 = 1.98733201817135256;
double convol::Q34 = 1.0;
double convol::MAX_SIZE_MASK_0 = 3.09023230616781;
double convol::MAX_SIZE_MASK_1 = 3.46087178201605;
double convol::MAX_SIZE_MASK_2 = 3.82922419517181;

int convol::MASK_SIZE(double MAX, double sigma) { return static_cast<int> (ceil(MAX * sigma)); }

convol::convol(int32_t cwidth, int32_t cheight) {
  width = cwidth;
  height = cheight;
}

double convol::normal(double x) {
  int sn;
  double R1, R2, y, y2, y3, y4, y5, y6, y7;
  double erf, erfc, z, z2, z3, z4;
  double phi;

  if (x < -UPPERLIMIT)
    return 0.0;
  if (x > UPPERLIMIT)
    return 1.0;

  y = x / SQRT2;
  if (y < 0) {
    y = -y;
    sn = -1;
  } else
    sn = 1;

  y2 = y * y;
  y4 = y2 * y2;
  y6 = y4 * y2;

  if (y < 0.46875) {
    R1 = P10 + P11 * y2 + P12 * y4 + P13 * y6;
    R2 = Q10 + Q11 * y2 + Q12 * y4 + Q13 * y6;
    erf = y * R1 / R2;
    if (sn == 1)
      phi = 0.5 + 0.5 * erf;
    else
      phi = 0.5 - 0.5 * erf;
  } else if (y < 4.0) {
    y3 = y2 * y;
    y5 = y4 * y;
    y7 = y6 * y;
    R1 = P20 + P21 * y + P22 * y2 + P23 * y3 + P24 * y4 + P25 * y5 + P26 * y6 +
         P27 * y7;
    R2 = Q20 + Q21 * y + Q22 * y2 + Q23 * y3 + Q24 * y4 + Q25 * y5 + Q26 * y6 +
         Q27 * y7;
    erfc = exp(-y2) * R1 / R2;
    if (sn == 1)
      phi = 1.0 - 0.5 * erfc;
    else
      phi = 0.5 * erfc;
  } else {
    z = y4;
    z2 = z * z;
    z3 = z2 * z;
    z4 = z2 * z2;
    R1 = P30 + P31 * z + P32 * z2 + P33 * z3 + P34 * z4;
    R2 = Q30 + Q31 * z + Q32 * z2 + Q33 * z3 + Q34 * z4;
    erfc = (exp(-y2) / y) * (1.0 / SQRTPI + R1 / (R2 * y2));
    if (sn == 1)
      phi = 1.0 - 0.5 * erfc;
    else
      phi = 0.5 * erfc;
  }

  return phi;
}

/** Integral of the Gaussian function */
double convol::phi0(double x, double sigma)

{
  return normal(x / sigma);
}

/** The Gaussian function */
double convol::phi1(double x, double sigma) {
  double t;

  t = x / sigma;
  return SQRT_2_PI_INV / sigma * exp(-0.5 * t * t);
}

/** First derivative of the Gaussian function */
double convol::phi2(double x, double sigma) {
  double t;

  t = x / sigma;
  return -x * SQRT_2_PI_INV / pow(sigma, 3.0) * exp(-0.5 * t * t);
}

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
void convol::compute_gauss_mask_0(double sigma, std::vector<double> &h) {
  int i, n;

  n = MASK_SIZE(MAX_SIZE_MASK_0, sigma); /* Error < 0.001 on each side */
  h.resize(2 * n + 1);

  for (i = -n + 1; i <= n - 1; i++)
    h[i + n] = phi0(-i + 0.5, sigma) - phi0(-i - 0.5, sigma);

  h[0] = 1.0 - phi0(n - 0.5, sigma);
  h[2 * n] = phi0(-n + 0.5, sigma);
}

/** First derivative of Gaussian smoothing mask */
void convol::compute_gauss_mask_1(double sigma, std::vector<double> &h) {
  int i, n;

  n = MASK_SIZE(MAX_SIZE_MASK_1, sigma); /* Error < 0.001 on each side */
  h.resize(2 * n + 1);

  for (i = -n + 1; i <= n - 1; i++)
    h[i + n] = phi1(-i + 0.5, sigma) - phi1(-i - 0.5, sigma);

  h[0] = -phi1(n - 0.5, sigma);
  h[2 * n] = phi1(-n + 0.5, sigma);
}

/** Second derivative of Gaussian smoothing mask */
void convol::compute_gauss_mask_2(double sigma, std::vector<double> &h) {
  int i, n;

  n = MASK_SIZE(MAX_SIZE_MASK_2, sigma); /* Error < 0.001 on each side */
  h.resize(2 * n + 1);

  for (i = -n + 1; i <= n - 1; i++)
    h[i + n] = phi2(-i + 0.5, sigma) - phi2(-i - 0.5, sigma);
  h[0] = -phi2(n - 0.5, sigma);
  h[2 * n] = phi2(-n + 0.5, sigma);
}

/** Convolve an image with the derivatives of a Gaussian smoothing kernel.
   Since all of the masks are separable, this is done in two steps in the
   function convolve_gauss.  Firstly, the rows of the image are convolved by
   an appropriate one-dimensional mask in convolve_rows_gauss, yielding an
   intermediate float-image h.  Then the columns of this image are convolved
   by another appropriate mask in convolve_cols_gauss to yield the final
   result k.  At the border of the image the gray values are mirrored. */

/** Convolve the rows of an image with the derivatives of a Gaussian. **/
void convol::convolve_rows_gauss(std::vector<double> &image,
                                 std::vector<double> &mask,
                                 std::vector<float> &h) {
  int32_t j, r, c, l;
  int n;

  n = (int)((mask.size() - 1) * 0.5);
  /* Inner region */
  for (r = n; r < height - n; r++) {
    for (c = 0; c < width; c++) {
      l = LCOR(r, c, width);
      float sum = 0.0;
      for (j = -n; j <= n; j++)
        sum += (double)(image[(l + j * width)]) * mask[j + n];
      h[l] = (float) sum;
    }
  }
  /* Border regions */
  for (r = 0; r < n; r++) {
    for (c = 0; c < width; c++) {
      l = LCOR(r, c, width);
      float sum = 0.0;
      for (j = -n; j <= n; j++)
        sum += (double)(image[LCOR(BR(r + j), c, width)]) * mask[j + n];
      h[l] = (float) sum;
    }
  }
  for (r = height - n; r < height; r++) {
    for (c = 0; c < width; c++) {
      l = LCOR(r, c, width);
      float sum = 0.0;
      for (j = -n; j <= n; j++)
        sum += (double)(image[LCOR(BR(r + j), c, width)]) * mask[j + n];
      h[l] = (float) sum;
    }
  }
}

/** Convolve the columns of an image with the derivatives of a Gaussian. */
void convol::convolve_cols_gauss(std::vector<float> &h,
                                 std::vector<double> &mask,
                                 std::vector<float> &k) {
  int32_t j, r, c, l;
  int n;

  n = (int)((mask.size() - 1) * 0.5);

  // Inner region
  for (r = 0; r < height; r++) {
    for (c = n; c < width - n; c++) {
      l = LCOR(r, c, width);
      float sum = 0.0;
      for (j = -n; j <= n; j++)
        sum += h[(l + j)] * mask[j + n];
      k[l] = (float)sum;
    }
  }
  // Border regions
  for (r = 0; r < height; r++) {
    for (c = 0; c < n; c++) {
      l = LCOR(r, c, width);
      float sum = 0.0;
      for (j = -n; j <= n; j++)
        sum += h[LCOR(r, BC(c + j), width)] * mask[j + n];
      k[l] = (float)sum;
    }
  }
  for (r = 0; r < height; r++) {
    for (c = width - n; c < width; c++) {
      l = LCOR(r, c, width);
      float sum = 0.0;
      for (j = -n; j <= n; j++)
        sum += h[LCOR(r, BC(c + j), width)] * mask[j + n];
      k[l] = (float)sum;
    }
  }
}

/** Convolve an image with a derivative of the Gaussian. */
void convol::convolve_gauss(std::vector<double> &image, double sigma,
                            int deriv_type, std::vector<float> &out) {
  std::vector<double> hr;
  std::vector<double> hc;

  switch (deriv_type) {
  case DERIV_R:
    compute_gauss_mask_1(sigma, hr);
    compute_gauss_mask_0(sigma, hc);
    break;
  case DERIV_C:
    compute_gauss_mask_0(sigma, hr);
    compute_gauss_mask_1(sigma, hc);
    break;
  case DERIV_RR:
    compute_gauss_mask_2(sigma, hr);
    compute_gauss_mask_0(sigma, hc);
    break;
  case DERIV_RC:
    compute_gauss_mask_1(sigma, hr);
    compute_gauss_mask_1(sigma, hc);
    break;
  case DERIV_CC:
    compute_gauss_mask_0(sigma, hr);
    compute_gauss_mask_2(sigma, hc);
    break;
  default: // just a stub
      assert(false);
    break;
  }

  // maskr = hr + nr;
  // maskc = hc + nc;
  std::vector<float> h(width * height, 0);
  convolve_rows_gauss(image, hr, h);
  convolve_cols_gauss(h, hc, out);
}

/** Convolve an image with a derivative of the Gaussian. */
void convol::convolve_gauss_opencv(std::vector<double>& image, double sigma,
    int deriv_type, std::vector<float>& out) {
    std::vector<double> hr;
    std::vector<double> hc;

    switch (deriv_type) {
    case DERIV_R:
        compute_gauss_mask_1(sigma, hr);
        compute_gauss_mask_0(sigma, hc);
        break;
    case DERIV_C:
        compute_gauss_mask_0(sigma, hr);
        compute_gauss_mask_1(sigma, hc);
        break;
    case DERIV_RR:
        compute_gauss_mask_2(sigma, hr);
        compute_gauss_mask_0(sigma, hc);
        break;
    case DERIV_RC:
        compute_gauss_mask_1(sigma, hr);
        compute_gauss_mask_1(sigma, hc);
        break;
    case DERIV_CC:
        compute_gauss_mask_0(sigma, hr);
        compute_gauss_mask_2(sigma, hc);
        break;
    default: // just a stub
        assert(false);
        break;
    }

    // @note: result is nearly identical to c version
    // if sepFilter2D is operates on double input and the output is transformed to float

    cv::Mat image_m(height, width, cv::DataType<double>::type, image.data());
    cv::Mat hr_m(hr.size(), 1, cv::DataType<double>::type, hr.data());
    cv::Mat hc_m(hc.size(), 1, cv::DataType<double>::type, hc.data());

    cv::Mat outm;
    cv::sepFilter2D(image_m, outm, -1, hc_m, hr_m, cv::Point(-1,-1), 0, CV_HAL_BORDER_REFLECT);

    std::vector<double> dimg(width * height);
    std::memcpy(dimg.data(), outm.ptr(0), outm.total() * sizeof(double));
    std::transform(dimg.begin(), dimg.end(), out.begin(), [](const double d) { return float(d); });

}

void convol::get_all_derivatives(std::vector<double>& image, double sigma,
    std::vector<std::vector<float>>& k) {

    static double zero_val = 0.0;
    k.resize(0);
    for (int i = 0; i < 5; i++)
        k.emplace_back(width * height, zero_val);

    timer::precise_stopwatch stopwatch;
    convolve_gauss_opencv(image, sigma, DERIV_R, k[0]);
    convolve_gauss_opencv(image, sigma, DERIV_C, k[1]);
    convolve_gauss_opencv(image, sigma, DERIV_RR, k[2]);
    convolve_gauss_opencv(image, sigma, DERIV_RC, k[3]);
    convolve_gauss_opencv(image, sigma, DERIV_CC, k[4]);

    if (debug()) {
        auto actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        std::cout << "Execution Time " << actual_wait_time << " microseconds" << std::endl;;
    }

}

/** Mirror the row coordinate at the borders of the image; height must be a
        defined variable in the calling function containing the image height. */
int32_t convol::BR(int32_t row) {
    return ((row) < 0 ? -(row) : (row) >= height ? height - (row)+height - 2 : (row));
}

/** Mirror the column coordinate at the borders of the image; width must be a
        defined variable in the calling function containing the image width. */
int32_t convol::BC(int32_t col) {
    return ((col) < 0 ? -(col) : (col) >= width ? width - (col)+width - 2 : (col));
}
