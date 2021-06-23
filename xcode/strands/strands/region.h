#pragma once
#include "chord.h"
#include <vector>

using namespace std;

namespace stegers
{
	/** Run-length encoded region of an image.  This type is returned by the
	threshold() function.  It provides the means to efficiently link line points
	into lines. */
	class region
	{
	public:
		region(const std::vector<int32_t>& image, uint32_t min_val, 
			int32_t image_width, int32_t image_m_height);

		int32_t num = 0; // number of chords
		std::vector<chord> rl; // array of chords

		static bool test();
	};

}


