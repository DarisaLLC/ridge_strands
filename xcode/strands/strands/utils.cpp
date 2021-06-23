
#include "utils.hpp"


std::vector<std::vector<double>> RectangularDoubleVector(int size1, int size2)
{
    static double zero = 0.0;
    std::vector<std::vector<double>> newVector(size1);
    for (int vector1 = 0; vector1 < size1; vector1++)
    {
        newVector[vector1] = std::vector<double>(size2, zero);
    }

    return newVector;
}


std::vector<std::vector<float>> RectangularFloatVector(int size1, int size2)
{
    std::vector<std::vector<float>> newVector(size1);
    for (int vector1 = 0; vector1 < size1; vector1++)
    {
        newVector[vector1] = std::vector<float>(size2);
    }

    return newVector;
}

/*
 * Utility function to draw a shape into a pelbuffer
 * optional scaling
 */
void DrawShape(cv::Mat& image, const char** shape, int32_t x, int32_t y, int32_t scale)
{
    int32_t width = x;
    int32_t m_height = y;
    for (auto s = shape; *s; s++, m_height += scale)
    {
        int32_t this_width = static_cast<int32_t>(strlen(*s)) * scale + x;
        if (this_width > width)
            width = this_width;
    }

    cv::Mat src(m_height, width, CV_8U);
    src.setTo(0);

    for (; *shape; shape++)
    {
        for (auto valPtr = *shape; *valPtr; valPtr++)
        {
            for (int32_t xx = 0; xx < scale; xx++)
            {
                for (int32_t yy = 0; yy < scale; yy++)
                {
                    int rr = y + yy;
                    int cc = static_cast<int>(x + (valPtr - *shape) * scale);
                    if (cc < 0 || cc >(width - 1) || rr < 0 || rr >(m_height - 1)) continue;

                    uint8_t* pel = src.ptr(rr,cc);
                    *pel = (*valPtr - '0');
                }
            }
        }
        y += scale;
    }

    image = src;
}
