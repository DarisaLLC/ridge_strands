#include <opencv2/opencv.hpp>
#include <iostream>
#include <boost/filesystem.hpp>
#include "derivatives.h"
#include "labelconnect.hpp"
#include "region.h"
#include "link.h"
#include "utils.hpp"
#include <iostream>
#include <iomanip>
#include "strands.h"
#include <cmath>
#include "version.h"
#include <exception> 
#include <time.h>

#include <chrono>

string NowToString()
{
	chrono::system_clock::time_point p = chrono::system_clock::now();
	time_t t = chrono::system_clock::to_time_t(p);
//	char str[26];
//	ctime_s(str, sizeof str, &t);
	return ctime(&t);
}

using namespace cv;
using namespace std;
using namespace boost;
namespace fs = boost::filesystem;
using namespace stegers;
using namespace timer;


void strands::run(){

	assert(labelConnect::test());
	assert(region::test());

	auto polarity = strands::line_polarity::light;
	float sigma = 1.5;
	float low_thr = 1.0;
	float high_thr = 6.0;

	std::string version_string(VERSION);
	version_string = NowToString() + "build " + version_string + "\nimage file: " + m_input_file.stem().string();
	auto filename = m_input_file.stem();
	filename = "histogram_output_" + filename.string() + ".csv";
	filename = m_output_dir / filename;
	
	Mat image = imread(m_input_file.string(), cv::IMREAD_GRAYSCALE);

	Mat display;
	string screen("display");
	if (show() || graphics()) {
		std::vector<Mat> triple{ image, image, image };
		merge(triple, display);
	}
	// Show our image inside a window.
	if (show())
		namedWindow(screen.c_str(), cv::WINDOW_GUI_EXPANDED);

	std::map<int32_t, int32_t> lm;

	try {


		// Convert to flat vector of doubles	
		// @todo have cv::Mat api
		cv::Mat dst;
		image.convertTo(dst, CV_64F);
		std::vector<double> dimg(size_t(image.rows * image.cols));
		assert(dst.isContinuous());
		std::memcpy(dimg.data(), dst.ptr(0), image.total() * sizeof(double));
		convol conv(image.cols, image.rows);
		conv.debug(debug());

		std::vector<std::vector<float>> out;
		conv.get_all_derivatives(dimg, sigma, out);
		std::vector<int32_t> ismax;
		std::vector<std::vector<double>> line_out;
		class link lnk(image.cols, image.rows, 1.5, strands::line_polarity::light);
		lnk.debug(debug());
		lnk.compute_line_points(out, ismax, line_out, low_thr, high_thr);
		const cv::Mat& esp = lnk.eSpace();
		lnk.compute_contours(ismax, line_out[0], line_out[5], line_out[6], line_out[1], line_out[2], line_out[3], line_out[5]);

		float length_total = 0;
		for (auto cc = 0; cc < lnk.results().size(); cc++) {
			const contour& ct = lnk.results()[cc];
			auto length = ct.compute_length();
			int ilength = (int)length;
			auto miter = lm.find(ilength);
			if (miter == lm.end())
				lm[ilength] = 0;
			auto current = lm[ilength];
			lm[ilength] = current + 1;
			length_total += length;
		}

		float avg_node_count = lnk.results().empty() ? 0 : length_total / float(lnk.results().size());
		std::string output = "AvgLength: " + std::to_string(avg_node_count);
		std::string output2 = "Count: " + std::to_string(lnk.results().size());

		if (debug()) {
			std::cout << output << std::endl;
			std::cout << output2 << std::endl;
		}

		string total_string = version_string + "\n" + output + "\n" + output2 + "\n";
		if (lnk.is_crowded())
			total_string = total_string + "is Crowded \n";

		save_csv(lm, filename.string(), total_string);




		if (show() || graphics()) {
			int factor = 8;
			for (auto cc = 0; cc < lnk.results().size(); cc++) {
				const contour& ct = lnk.results()[cc];
				std::vector<Point> pts;
				std::vector<float>::const_iterator rowItr = ct.row.begin();
				std::vector<float>::const_iterator colItr = ct.col.begin();
				for (auto pp = 0; pp < ct.row.size(); pp++, rowItr++, colItr++) {
					pts.emplace_back(static_cast<int>(*colItr * factor), static_cast<int>(*rowItr * factor));
				}
				polylines(display, pts, false, Scalar(0, 0, 255), 2, LINE_AA, int(std::log2(factor)));
			}


			if (lnk.is_crowded())
				putText(display, "Is Crowded ", Point(1000, 900), FONT_HERSHEY_SIMPLEX, 3.0, Scalar(0, 255, 0), 2);


			std::string output = " Average Length " + std::to_string(avg_node_count);
			putText(display, output.c_str(), Point(1000, 1000), FONT_HERSHEY_SIMPLEX, 3.0, Scalar(0, 255, 0), 2);
			std::string output2 = " Count " + std::to_string(lnk.results().size());
			putText(display, output2.c_str(), Point(1000, 1200), FONT_HERSHEY_SIMPLEX, 3.0, Scalar(0, 255, 0), 2);

			int start = 100;
			auto it = lm.begin();
			while (it != lm.end()) {
				std::string hout = "histogram [" + std::to_string(it->first) + "] = " + std::to_string(it->second);
				putText(display, hout.c_str(), Point(500, start), FONT_HERSHEY_SIMPLEX, 1.0, Scalar(255, 0, 0), 2);
				start += 30;
				it++;
			}

			if (show()) {
				imshow(screen.c_str(), display);
				// Wait for any keystroke in the window
				waitKey(0);
			}

			if (graphics()) {
				auto filename = m_input_file.filename();
				filename = "output_" + filename.string();
				filename = m_output_dir / filename;
				cv::imwrite(filename.string(), display);
			}

		}

	}

	catch (std::exception& e) {
		version_string = version_string + " error " + e.what();
		lm.clear();
		save_csv(lm, filename.string(), version_string);
	}

}

