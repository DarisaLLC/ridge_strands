#include "chord.h"

namespace stegers
{

	chord::chord()
	{
			  r = 0;
			  cb = 0;
			  ce = 0;
	}

	chord::chord(short nr, short ncb, short nce)
	{
		  r = nr;
		  cb = ncb;
		  ce = nce;
	}
}
