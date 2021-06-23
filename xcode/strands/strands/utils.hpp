#pragma once
#pragma once
#include <vector>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <boost/filesystem.hpp>

#include <chrono>
#include <atomic>

using namespace boost;
namespace bfs = boost::filesystem;

// https://codereview.stackexchange.com/a/196253

namespace timer
{
    template <typename Clock = std::chrono::high_resolution_clock>
    class stopwatch
    {
        const typename Clock::time_point start_point;
    public:
        stopwatch() :
            start_point(Clock::now())
        {}

        typedef typename Clock::duration::rep Rep;

        template <typename Rep = typename Clock::duration::rep, typename Units = typename Clock::duration>
        Rep elapsed_time() const
        {
            std::atomic_thread_fence(std::memory_order_relaxed);
            auto counted_time = std::chrono::duration_cast<Units>(Clock::now() - start_point).count();
            std::atomic_thread_fence(std::memory_order_relaxed);
            return static_cast<Rep>(counted_time);
        }
    };

    using precise_stopwatch = stopwatch<>;
    using system_stopwatch = stopwatch<std::chrono::system_clock>;
    using monotonic_stopwatch = stopwatch<std::chrono::steady_clock>;
}


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#define LCOR(row,col,width) (row)*(width) + (col)

std::vector<std::vector<double>> RectangularDoubleVector(int size1, int size2);
std::vector<std::vector<float>> RectangularFloatVector(int size1, int size2);

/*
 * Utility function to draw a shape into a pelbuffer
 * optional scaling
 */
void DrawShape(cv::Mat& image, const char** shape, int32_t x = 0, int32_t y = 0, int32_t scale = 1);


/*
 * RAII (Resource Allocation is Initialization) Exception safe handling of openning and closing of files 
 */

struct ofStreamDeleter {
    void operator()(std::ofstream* p) const {
        if (p != nullptr) {
            p->close();
            delete p;
        }
    }
};

static std::shared_ptr<std::ofstream> make_shared_ofstream(std::ofstream* ifstream_ptr) {
    return std::shared_ptr<std::ofstream>(ifstream_ptr, ofStreamDeleter());
}

static std::shared_ptr<std::ofstream> make_shared_ofstream(const std::string& filename) {
    return make_shared_ofstream(new std::ofstream(filename, std::ofstream::out));
}


template<class T, template<typename ELEM, typename ALLOC=std::allocator<ELEM>> class CONT=std::vector>
static void save_csv(const CONT<T>& data, const bfs::path& output_file) {
    static std::string delim(",");
    auto papa = output_file.parent_path();
    if (bfs::exists(papa)) {
        std::shared_ptr<std::ofstream> myfile = make_shared_ofstream(output_file.string());
        auto cnt = 0;
        auto size = data.size() - 1;
        for (const T& dd : data) {
            *myfile << dd;
            if (cnt++ < size)
                *myfile << delim;
        }
        *myfile << std::endl;
    }
}

template <typename Key, typename T>
static void save_csv(const std::map<Key, T>& m, const bfs::path& output_file, const std::string& prefix) {
    static std::string delim(",");
    auto papa = output_file.parent_path();
    if (bfs::exists(papa)) {
        std::shared_ptr<std::ofstream> myfile = make_shared_ofstream(output_file.string());
        *myfile << prefix << std::endl;

        for (typename std::map<Key, T>::const_iterator it = m.begin();
            it != m.end(); it++) {
            *myfile << it->first << delim << it->second << std::endl;
        }
        *myfile << std::endl;
    }
}

