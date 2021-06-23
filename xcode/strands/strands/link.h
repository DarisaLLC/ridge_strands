#pragma once
#include <assert.h>
#include <ostream>
#include <vector>
#include "derivatives.h"
#include <boost/math/constants/constants.hpp>
#include <opencv2/opencv.hpp>
#include "contour.h"
#include "strands.h"

using namespace std;
using namespace cv;

const double pi = boost::math::constants::pi<double>();

class crossRef {
    public:
    /*
    Storage the Crossref variables, it is the Correction.java code
    This data structure facilitates the quick search for the next possible starting point of a line.An array of crossrefs will be accumulatedand
    sorted according to its value.xand y are the coordinates of a point in the image.When this point has been processed it will be marked as done.
    */
    crossRef(int32_t x = 0, int32_t y = 0, double value = 0.0, bool done = false) : m_x(x), m_y(y), m_val(value), m_done(done) {}

    // Accessors
    const int32_t& x() const { return m_x; }
    const int32_t& y() const { return m_y; }
    const double& value() const { return m_val;  }
    bool done() const { return m_done;  }

    void setDone() const { m_done = true; }
    void setUnDone() const { m_done = false; }

    int32_t compareTo(crossRef& other) {
        int32_t rt = (m_val > other.value()) ? -1 : (m_val < other.value()) ? 1 : 0;
        return rt;
    }

    bool operator==(const crossRef& other) const {
        return this->m_val == other.value();
    }


    bool operator<(const crossRef& other) const {
        return this->m_val < other.value();
    }


    bool operator>(const crossRef& other) const {
        return this->m_val > other.value();
    }

    friend ostream& operator<< (ostream& ous, const crossRef& dis)
    {
        //    return  + str(self.x) +  + str(self.y) +  + str(self.value) +  + str(self.done)
        ous << "x: " << dis.x() << "\ty: " << dis.y() << "\tvalue: " << dis.value() << "\tdone: " << std::boolalpha << dis.done();
        return ous;
    }


    private:
        int32_t m_x, m_y;
        double m_val;
		mutable bool m_done;

        };





	/** Data structure to store three doubles: x,y and t (distance along line) */
class doublepoint
{

public:
    double cx = 0;
    double cy = 0;
    double t = 0;

    doublepoint() {}

    doublepoint(double ccx, double ccy, double ct) {
            cx = ccx;
            cy = ccy;
            t = ct;
    }

    doublepoint(double ccx, double ccy) {
		cx = ccx;
		cy = ccy;
		t = 0;
    }
};

/** Offsets to a specific location in the image.  An array of this type is
returned by the modified Bresenham algorithm in width.c.  It is also used
in link.c to hold an array of pixel locations to check for appropriate
neighbors. */
class offset
{
public:
    int32_t x;
    int32_t y;

    offset() { x = 0; y = 0; }

    offset(int32_t nx, int32_t ny) {
        x = nx;
        y = ny;
    }

};


/** This data structure is used to accumulate junction information.  It is
needed to split lines at junction points. */
class junction
{
    /** Index of line that is already processed */
public:
    int32_t cont1 = 0;
    /** Index of line that runs into cont1 */
    int32_t cont2 = 0;
    /** Index of the junction point in cont1 */
    int32_t pos = 0;
    /** y-(row-)coordinate of the junction point (corrected for ImageJ)*/
    float x = 0;
    /** x-(col-)coordinate of the junction point (corrected for ImageJ)*/
    float y = 0;

    junction() {}

    junction(int32_t ncont1, int32_t ncont2, int32_t npos, float nx, float ny) {
        cont1 = ncont1;
        cont2 = ncont2;
        pos = npos;
        x = nx;
        y = ny;
    }


    bool operator==(const junction& other) const {
        return this->pos == other.pos;
    }


    bool operator<(const junction& other) const {
        return this->pos < other.pos;
    }


    bool operator>(const junction& other) const {
        return this->pos > other.pos;
    }

    /** This function compares two junctions according to their first line indexes,
     and, if needed, by the position of the junction within the line.  It is
     called by qsort. */
    static bool compare(const junction& thisJunction, const junction& otherjunction) {
        if (thisJunction.cont1 == otherjunction.cont1)
        {
            return thisJunction.pos >= otherjunction.pos;
        }
        else
        {
            return thisJunction.cont1 > otherjunction.cont1;
        }
    }


    friend ostream& operator<< (ostream& ous, const junction& dis)
    {
        ous << "cont1: " << dis.cont1 << "\tcont2: " << dis.cont2 << "\tpos: " << dis.pos << "\txy: " << dis.x << "," << dis.y;
        return ous;
    }
};

class link {
public:
    link(int32_t width, int32_t height, float sigma = 1.0,
        strands::line_polarity lp = strands::line_polarity::dark, int min_node_count = 5,
        bool extend_lines = true, bool split_lines = false)
        : m_width(width), m_height(height), m_extend_lines(extend_lines), m_sigma(sigma),
        m_polarity(lp), m_split_lines(split_lines), nMinNumberOfNodes(min_node_count){
        m_total_cnt = m_width * m_height;

        // Parameters @todo move into params class
        MAX_ANGLE_DIFFERENCE = pi / 6.0;
        MAX_LINE_EXTENSION = m_sigma * 2.5;
        MAX_LINE_WIDTH = (2.5 * m_sigma);
        LINE_WIDTH_COMPENSATION = 1.05;
        MIN_LINE_WIDTH = 0.1;
        m_pixel_boundary = 0.6;
        m_spread_fraction = 0.1;
        m_espace = cv::Mat (m_height, m_width, CV_32F);
        m_max_line_pts = (m_width * m_height) / (m_sigma *  9);
        m_max_line_pts = m_spread_fraction * m_max_line_pts;
        m_crowded = false;
    }

    inline const int32_t& width() const { return m_width; }
    inline const int32_t& height() const { return m_height; }
    inline const int32_t& total() const { return m_total_cnt; }

    const strands::line_polarity& polarity() const { return m_polarity;  }
    bool extend_lines() const { return m_extend_lines; }
    bool split_lines() const { return m_split_lines; }
    float sigma() const { return m_sigma; }
    int minNodeCount () const{ return nMinNumberOfNodes;  }
    bool is_crowded() const { return m_crowded;  }


    const std::vector<contour>& results() const { return m_results; }
    const std::vector<junction>& junctions() const { return m_junctions; }

    bool debug() const { return m_debug; }
    void debug(bool state) const { m_debug = state;  }


    void compute_line_points(const std::vector<std::vector<float>>& derivs, std::vector<int32_t>& ismax,
        std::vector<std::vector<double>>& enpp, double low, double high);

    void compute_contours(std::vector<int32_t> ismax,
        const std::vector<double>& eigval,
        const std::vector<double>& posx, const std::vector<double>& posy,
        const std::vector<double>& normx, const std::vector<double>& normy,
        const std::vector<double>& gradx, const std::vector<double>& grady);

    const cv::Mat& eSpace() const { return m_espace;  }

private:

    std::vector<contour> m_results;
    std::vector<junction> m_junctions;
    mutable cv::Mat m_espace;

    // orintation from norms. 
    double get_angle(const std::vector<double>& nrmx, const std::vector<double>& nrmy, int32_t position) {

        auto nx = -nrmy[position];
        auto ny = nrmx[position];
        auto alpha = std::atan2(ny, nx);
        if (alpha < 0.0)
            alpha += 2.0 * pi;
        if (alpha >= pi)
            alpha -= pi;
        return alpha;
    }


    int32_t BR(int32_t row) {
        return ((row) < 0 ? -(row) : (row) >= m_height ? m_height - (row)+m_height - 2 : (row));
    }

    /** Mirror the column coordinate at the borders of the image; width must be a
            defined variable in the calling function containing the image width. */
    int32_t BC(int32_t col) {
        return ((col) < 0 ? -(col) : (col) >= m_width ? m_width - (col)+m_width - 2 : (col));
    }


    void getismax(double low, double high, const vector<double>& eigval, const vector<double>& posx, const vector<double>& posy, std::vector<uint32_t>& iismax);
    void bresenham(double nx, double ny, double px, double py, double length, std::vector<offset>& out);

    double interpolate_response(const std::vector<double>& resp, int32_t x, int32_t y, double px, double py);

/** Calculate the closest point to (px,py) on the line (lx,ly) + t*(dx,dy)
 *  and return the result in (cx,cy), plus the parameter in t. */
	void closest_point(double lx, double ly, double dx, double dy, double px, double py, doublepoint& out);

/**
 * Interpolate the gradient of the gradient images gradx and grady with width
 * width at the point (px,py) using linear interpolation, and return the
 * result in (gx,gy) (doublepoint format)
 * */
	void interpolate_gradient(const std::vector<double>& gradx, const std::vector<double>& grady, double px, double py, doublepoint& out);


    mutable int32_t m_width, m_height, m_total_cnt;
    mutable strands::line_polarity m_polarity;
    mutable bool m_extend_lines, m_split_lines;
    mutable bool m_debug;

    float m_sigma;
    double m_pixel_boundary;
    float m_spread_fraction;
    size_t m_max_line_pts;
    bool m_crowded;
    double MAX_ANGLE_DIFFERENCE;
    double MAX_LINE_EXTENSION;
    double MAX_LINE_WIDTH;
    double LINE_WIDTH_COMPENSATION;
    double MIN_LINE_WIDTH;
    double MAX_CONTRAST;

    int  nMinNumberOfNodes;
};
