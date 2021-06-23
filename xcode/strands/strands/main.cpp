
#include <fstream>
#include <iostream>
#include <iomanip>  //  setw
#include <cmath>  //  floor
#include <climits>  //  UINT_MAX

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>   // includes all needed Boost.Filesystem declarations
#include "strands.h"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

using namespace std;
namespace bf = boost::filesystem;


int main(int argc, char* argv[]) {
	po::options_description desc("Usage");
	desc.add_options()
	("help", "this help message")
	("input,I", po::value< string >(), "image file path ")
	("output,O", po::value< string >(), "output directory path")
	("show,S", po::value< bool >(), "Display Results")
	("graphics,G", po::value< bool >(), "Write Results Display")
	("debug,D", po::value< bool >(), "output debugging info")
	;
	
	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	}
	catch (std::exception& e) {
		cout << "Invalid parameter: " << e.what() << endl;
		return 1;
	}
	
	if (vm.count("help")) {
			// cout << "Melter version " << STRANDS_VERSION_MAJOR << "." << STRANDS_VERSION_MINOR << endl;
		cout << desc << endl;
		return 1;
	}
	
	bf::path input_path;
	if (vm.count("input")) {
		try {
			input_path = vm["input"].as< string >();
		}
		catch (std::exception& e) {
			cout << "Invalid parameter: " << e.what() << endl;
			return 1;
		}
	}
	else {
		cout << "No input file provided." << endl;
		return 1;
	}
	
	bf::path output_path;
	if (vm.count("output")) {
		try {
			output_path = vm["output"].as< string >();
		}
		catch (std::exception& e) {
			cout << "Invalid parameter: " << e.what() << endl;
			return 1;
		}
	}
	
	bool show_on = false;
	if (vm.count("show")) {
		try {
			show_on = vm["show"].as< bool >();
		}
		catch (std::exception& e) {
			cout << "Invalid parameter: " << e.what() << endl;
			return 1;
		}
	}
	
	bool output_graphics_on = false;
	if (vm.count("graphics")) {
		try {
			output_graphics_on = vm["graphics"].as< bool >();
		}
		catch (std::exception& e) {
			cout << "Invalid parameter: " << e.what() << endl;
			return 1;
		}
	}
	
	bool debug_on = false;
	if (vm.count("debug")) {
		try {
			debug_on = vm["debug"].as< bool >();
		}
		catch (std::exception& e) {
			cout << "Invalid parameter: " << e.what() << endl;
			return 1;
		}
	}
	
	if (debug_on) {
		cout << "Input file   : " << setw(32) << right << input_path << endl;
		cout << "Output file  : " << setw(32) << right << output_path << endl;
	}
	
	if (!bf::exists(input_path) || !bf::is_regular_file(input_path)) {
		cout << input_path << " does not exists or is a regular file." << endl;
		return 1;
	}
	
	if (!bf::exists(output_path) || !bf::is_directory(output_path)) {
		cout << output_path << " does not exists or is a directory." << endl;
		return 1;
	}
	
		// All good create a strands app
	strands stranded(input_path, output_path);
	stranded.show(show_on);
	stranded.graphics(output_graphics_on);
	stranded.debug(debug_on);
	
	stranded.run();
	
	return 0;
}
