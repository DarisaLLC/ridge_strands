#include "region.h"
#include "utils.hpp"

namespace stegers
{
	region::region(const std::vector<int32_t>& image, uint32_t min_val,
        int32_t image_width, int32_t image_m_height){
		rl.resize(0);

        long   grey;
        long   r, c, l, count;
        bool   inside;

        inside = false;
        count = 0;
        rl.emplace_back();

        for (r = 0; r < image_m_height; r++) {
            for (c = 0; c < image_width; c++) {
                l = LCOR(r, c, image_width);
                grey = image[l];
                if (grey >= min_val) {
                    if (!inside) {
                        inside = true;
                        rl[count].r = (int16_t)r;
                        rl[count].cb = (int16_t)c;
                    }
                }
                else {
                    if (inside) {
                        inside = false;
                        rl[count].ce = (int16_t)(c - 1);
                        count++;
                        rl.emplace_back();
                    }
                }
            }
            if (inside) {
                inside = false;
                rl[count].ce = (int16_t)(image_width - 1);
                count++;
                rl.emplace_back();
            }
        }
        this->num = count;
	}




    bool region::test() {
        const char* frame[] =
        {
            "00100100",
            "00110100",
            "00011000",
            "01000100",
            "01000000",
            0 };
        const char* gold[] =
        {
            "00100100",
            "00110100",
            "00011000",
            "02000100",
            "02000000",
            0 };

        cv::Mat pels(5, 8, CV_8U);
        DrawShape(pels, frame);

        std::vector<int32_t> lpels(5 * 8);
        for (auto row = 0; row < pels.rows; row++)
            for (auto col = 0; col < pels.cols; col++) {
                auto l = LCOR(row, col, pels.cols);
                lpels[l] = int32_t(pels.at<uint8_t>(row, col));
            }

        region rg(lpels, 1, pels.cols, pels.rows);
        bool check = rg.rl.size() == 9;
        if (!check) return check;

        std::vector<chord> golds;
        golds.emplace_back(0, 2, 2);
        golds.emplace_back(0, 5, 5);
        golds.emplace_back(1, 2, 3);
        golds.emplace_back(1, 5, 5);
        golds.emplace_back(2, 3, 4);
        golds.emplace_back(3, 1, 1);
        golds.emplace_back(3, 5, 5);
        golds.emplace_back(4, 1, 1);
        golds.emplace_back(0, 0, 0);

        for (int i = 0; i < 9; i++)
            assert(rg.rl[i] == golds[i]);

        return true;
    }
}