
#include "labelconnect.hpp"
#include "utils.hpp"

//----------------------------------------------------------------------------
// Connected component labeling.
//
// The 1D connected components of each row are labeled first.  The labels
// of row-adjacent components are merged by using an associative memory
// scheme.  The associative memory is stored as an array that initially
// represents the identity permutation M = {0,1,2,...}.  Each association of
// i and j is represented by applying the transposition (i,j) to the array.
// That is, M[i] and M[j] are swapped.  After all associations have been
// applied during the merge step, the array has cycles that represent the
// connected components.  A relabeling is done to give a consecutive set of
// positive component labels.
//
// For example:
//
//     Image       Rows labeled
//     00100100    00100200
//     00110100    00330400
//     00011000    00055000
//     01000100    06000700
//     01000000    08000000
//
// Initially, the associated memory is M[i] = i for 0 <= i <= 8.  Background
// label is 0 and never changes.
// 1. Pass through second row.
//    a. 3 is associated with 1, M = {0,3,2,1,4,5,6,7,8}
//    b. 4 is associated with 2, M = {0,3,4,1,2,5,6,7,8}
// 2. Pass through third row.
//    a. 5 is associated with 3, M = {0,3,4,5,2,1,6,7,8}.  Note that
//       M[5] = 5 and M[3] = 1 have been swapped, not a problem since
//       1 and 3 are already "equivalent labels".
//    b. 5 is associated with 4, M = {0,3,4,5,1,2,6,7,8}
// 3. Pass through fourth row.
//    a. 7 is associated with 5, M = {0,3,4,5,1,7,6,2,8}
// 4. Pass through fifth row.
//    a. 8 is associated with 6, M = {0,3,4,5,1,7,8,2,6}
//
// The M array has two cycles.  Cycle 1 is
//   M[1] = 3, M[3] = 5, M[5] = 7, M[7] = 2, M[2] = 4, M[4] = 1
// and cycle 2 is
//   M[6] = 8, M[8] = 6
// The image is relabeled by replacing each pixel value with the index of
// the cycle containing it.  The result for this example is
//     00100100
//     00110100
//     00011000
//     02000100
//     02000000
//----------------------------------------------------------------------------

#include <vector>

class u8s32LUT
{
public:
    u8s32LUT ()
    {
        lut_.resize(256, 1);
        lut_[0] = 0;
    }

    void map (uint8_t* src, int32_t* loc)
    {
        *loc = lut_[*src];
    }
    
    const std::vector<int32_t>& lut () const { return lut_; }
    
private:
    std::vector<int32_t> lut_;
};

#define SINGLETONPTR(Typename)                     \
  static Typename* instance() {                 \
    static Typename e;                          \
    return &e;                                  \
  }

SINGLETONPTR(u8s32LUT);



labelConnect::labelConnect(const cv::Mat& src)
{
    src_shallow_ = src;
}

const typename labelConnect::label_bbox_map_t& labelConnect::label_bbox_map () const { return label2roi_; }

const cv::Mat& labelConnect::label () const { return label_; }

uint32_t labelConnect::regions () const { return regions_; }

bool labelConnect::run () { get_components(); return true; }


void labelConnect::get_components ()
{
    // Create a temporary copy of image to store intermediate information
    // during component labeling.  The original image is embedded in an
    // image with two more rows and two more columns so that the image
    // boundary pixels are properly handled.
    // To prevent duplicate clears, when src is 32 bit we are assuming that it is 2 pixels
    // larger already. This is a hack to be cleanly handled
    
   
    const std::vector<int32_t>& lut = instance()->lut();
    
    buffer_ = cv::Mat(src_shallow_.rows + 2, src_shallow_.cols + 2, CV_32S);
    buffer_.setTo(0);
    label_ = buffer_(cv::Range(1, buffer_.rows-1), cv::Range(1, buffer_.cols-1));
    
    for (auto iY = 0; iY < src_shallow_.rows; ++iY)
        {
            const uint8_t *row = src_shallow_.ptr (iY);
            int32_t *brow = (int32_t *) label_.ptr (iY);
            for (auto iX = 0; iX < src_shallow_.cols; ++iX)
            {
                *brow++ = lut[*row++];
            }
        }

    // label connected components in 1D array
    icomponent_ = 0;
    
    for (auto iY = 1; iY < buffer_.rows - 1; ++iY)
    {
        // Next row, second column
        int32_t *brow = (int32_t *) buffer_.ptr (iY);
        ++brow;
        
        // Label runs with an incrementing component number
        for (auto iX = 1; iX < buffer_.cols - 1; ++iX)
        {
            if (brow[iX])
            {
                ++icomponent_;
                while ((iX < (buffer_.cols - 1)) && brow[iX])
                {
                    brow[iX] = icomponent_;
                    ++iX;
                }
            }
        }
    }
    
    // If there are no non zero pels
    // Labled component image is not assigned a real one
    if ( icomponent_ == 0 )
    {
        regions_ = 0;
        return;
    }
    
    // associative memory for merging
    amm_.resize(icomponent_+1);
    for (int32_t i = 0; i < icomponent_ + 1; ++i)
        amm_[i] = i;
    
    // Merge equivalent components.  Pixel (x,y) has previous neighbors
    // (x-1,y-1), (x,y-1), (x+1,y-1), and (x-1,y) [4 of 8 pixels visited
    // before (x,y) is visited, get component labels from them].
    
    const int32_t rup = buffer_.cols;
    
    for (auto iY = 1; iY < buffer_.rows-1; ++iY)
    {
        int32_t *brow = (int32_t *) buffer_.ptr (iY);
        ++brow;
        
        for (auto iX = 1; iX < buffer_.cols-1; ++iX)
        {
            const int32_t iValue = *brow++; // at +1 after this
            if ( iValue > 0 )
            {
                add_to_associative(iValue, brow [-rup - 2]); // brow[-rup - 2]
                add_to_associative(iValue, brow [-rup - 1]); // brow [-rup - 1]
                add_to_associative(iValue, brow [-rup]); // brow [-rup]
                add_to_associative(iValue, brow [rup - 2]); // brow [rup - 2]
            }
        }
    }
    
    // replace each cycle of equivalent labels by a single label
    regions_ = 0;
    for (auto i = 1; i <= icomponent_; ++i)
    {
        if ( i <= amm_[i] )
        {
            regions_++;
            uint32_t iCurrent = i;
            while ( amm_[iCurrent] != i )
            {
                const uint32_t iNext = amm_[iCurrent];
                amm_[iCurrent] = regions_;
                iCurrent = iNext;
            }
            amm_[iCurrent] = regions_;
        }
    }
    
    // pack a relabeled image in smaller size output
    // create a rcWindow at 1,1 and label according to the associative memory TBD
    
    for (auto iY = 0; iY < label_.rows; ++iY)
    {
        int32_t* brow = (int32_t *) label_.ptr (iY);
        
        for (auto iX = 0; iX < label_.cols; ++iX, ++brow)
        {
            *brow = amm_[*brow]; // *brow++ here appears to be a compiler bug
        }
    }
}


void labelConnect::get_rects ()
{
    int nLabels = regions_;
    
    // Regardless of label is a window or not, lower right should be initialized to zero
    const int m_height = static_cast<int>( label_.rows );
    const int width = static_cast<int>( label_.cols );
    iPair zeropair (0,0);
    std::vector<iPair> lrp (nLabels+1,zeropair);
    std::vector<iPair> ulp(nLabels + 1, iPair(label_.cols, label_.rows));
    
    /*
     * The following computes the min containing rectangle for
     * the connected region. That is for the integer top left position
     * of the most top left (lowest x and lowest y) pixel position in the runs and
     * the most bottom right (highest x and highest y) pixel position in
     * the runs.
     * The integer bounding box is at upper left x and y and has width and m_height
     * of lower - upper + 1.
     */
    
    for (auto j = 0; j < m_height; ++j)
    {
        int32_t* pSrc = (int32_t*) label_.ptr(j);
        for (auto i=0; i < width; ++i, ++pSrc)
        {
            const int32_t l = *pSrc;
            // label 0 is background. Currently we do not do any thing with that
            // Note: first rect is at index 1 to avoid computing l-1 inside the inner loop
            if (l)
            {
                // Update min enclosing rectangle
                if (i < ulp[l].x())
                    ulp[l].x() = i;
                if (j < ulp[l].y())
                    ulp[l].y() = j;
                
                if (i > lrp[l].x())
                    lrp[l].x() = i;
                if (j > lrp[l].y())
                    lrp[l].y() = j;
            }
        }
    }
  
    
    // Construct map of rects
    for (auto l = 1; l <= nLabels; ++l )
    {
        /* The integer bounding box is at upper left x and y and has width and m_height
         * of lower - upper + 1.
         * We use the rectangle ctor for passing the position and width and m_height
         */
        iRect r(ulp[l].x(), ulp[l].y(), lrp[l].x() - ulp[l].x() + 1, lrp[l].y() - ulp[l].y() + 1);

        label2roi_[l] = r;
    }
    
}


void labelConnect::add_to_associative (labelConnect::bin_type i0, labelConnect::bin_type i1)
{
    
    if ( i1 == 0 )
        return;
    
    int32_t iSearch = i1;
    do
    {
        iSearch = amm_[iSearch];
    }
    while ( iSearch != i1 && iSearch != i0 );
    
    if ( iSearch == i1 )
    {
        uint32_t iSave = amm_[i0];
        amm_[i0] = amm_[i1];
        amm_[i1] = iSave;
    }
}





bool labelConnect::test() {
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

    cv::Mat label(5, 8, CV_8U);
    DrawShape(label, gold);

    labelConnect lc(pels);
    lc.run();
    bool check = label.size() == lc.label().size();
    if (!check) return check;


    for (auto j = 0; j < label.rows; j++)
        for (auto i = 0; i < label.cols; i++)
        {
            auto this_check = (label.at<uint8_t>(j,i) == lc.label().at<int32_t>(j,i));
            if (!this_check) return this_check;
        }



    simpleRLE srlep8(label);
    {
        std::vector<simpleRLE::run_t> rles;
        int rle_count = srlep8.rle_row((uint8_t*)label.ptr(0), label.cols, rles);
        check &= rle_count == 5;
        if (!check) return check;
        check &= rles[0].first == 0;
        if (!check) return check;
        check &= rles[0].second == 2;
        if (!check) return check;
        check &= rles[1].first == 1;
        if (!check) return check;
        check &= rles[1].second == 1;
        if (!check) return check;
        check &= rles[2].first == 0;
        if (!check) return check;
        check &= rles[2].second == 2;
        if (!check) return check;
        check &= rles[3].first == 1;
        if (!check) return check;
        check &= rles[3].second == 1;
        if (!check) return check;
        check &= rles[4].first == 0;
        if (!check) return check;
        check &= rles[4].second == 2;
        if (!check) return check;
    }
    return check;
}






