

#ifndef labelconnect_hpp
#define labelconnect_hpp

#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <map>
#include "rectangle.hpp"

class labelConnect
{
public:
    typedef uint32_t bin_type;
    typedef std::map<bin_type,iRect> label_bbox_map_t;
    labelConnect (const cv::Mat& src);
    bool run ();
    const cv::Mat& label () const;
    uint32_t regions () const;
    const label_bbox_map_t& label_bbox_map () const;
    static bool test();
private:
    
    void get_components ();
    void get_rects ();
    cv::Mat src_shallow_;
    cv::Mat buffer_;
    cv::Mat label_;
    std::vector<iRect> rois_;
    std::vector<bin_type> amm_;
    bin_type regions_;
    bin_type icomponent_;
    label_bbox_map_t label2roi_;
    void add_to_associative (bin_type, bin_type);
};

//
//  Run Length Encoding of Mask Areas.
//
//
//
//
//


class simpleRLE
{
public:
    typedef uint8_t* pixel_ptr_t;
    typedef uint8_t pixel_t;
    typedef std::pair<pixel_t, uint16_t> run_t;

    // Rle all
    simpleRLE(const cv::Mat& src)
    {
        rat_.resize (0);
        for (auto row = 0; row < src.rows; row++)
        {
            std::vector<run_t> rles;
            rle_row ((uint8_t*) src.ptr(row), src.cols, rles);
            rat_.emplace_back(rles);
        }
    }
    
    // run_vector_t& pointToRow(int32_t y) const {return rat_[y];}
    
    int rle_row (pixel_ptr_t input, int length, std::vector<std::pair<pixel_t, uint16_t> >& output)
    {
        static const uint32_t max_run = numeric_limits<uint16_t>::max();
        
        int count = 0,index;
        pixel_t pixel;
        output.resize(0);
        int out = 0;
        
        while (count < length)
        {
            index = count;
            pixel = input[index++];
            while (index < length && index - count < max_run && input[index] == pixel)
                index++;
            output.emplace_back (pixel, (uint16_t)(index - count));
            count=index;
            out++;
        } /* while */
        return(out);
    }
    
    
    
private:
    mutable std::vector<std::vector<run_t > > rat_;

    
    };
    

    
#endif /* labelconnect_hpp */
