#ifndef BUILD_STRANDS_HPP
#define BUILD_STRANDS_HPP

#include <boost/filesystem.hpp>
using namespace boost;
namespace fs = boost::filesystem;
#include <string>
using namespace std;


#define STRANDS_VERSION_MAJOR 0
#define STRANDS_VERSION_MINOR 1


/*!
     The strands class is the main class for this program.  The main driver
     (in main.cpp) creates an instance of this class as the first task.  It is
     also the last class destroyed before the program exits.
*/
class strands {
public:


    enum line_polarity {
        light,
        dark
    };

    strands(const fs::path& image_file, const fs::path& output_dir, const std::string& result_file_name = "") :
        m_input_file(image_file), m_output_dir(output_dir), debug_flag(false), output_graphics_flag(false), show_flag(false) {
        m_result_filename = (result_file_name == "") ? m_input_file.stem().generic_string() : result_file_name;
    }

    bool show() const { return show_flag;  }
    void show(bool state) const { show_flag = state;  }
    bool graphics() const { return output_graphics_flag; }
    void graphics(bool state) const { output_graphics_flag = state; }
    bool debug() const { return debug_flag; }
    void debug(bool state) const { debug_flag = state; }
    void run();


private:
    //!  Set to true if debug output is required; false otherwise
    mutable bool debug_flag;
    //!  Set to true if verbose output is required; false otherwise
    mutable bool output_graphics_flag;
    //!  Set to true if graphical output is required; false otherwise
    mutable bool show_flag;

    fs::path m_input_file;
    fs::path m_output_dir;
    std::string m_result_filename;

    //!  Number of rows
    unsigned int M;
    //!  Number of columns
    unsigned int N;
};

#endif
