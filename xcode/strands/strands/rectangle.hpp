#ifndef __rcRECTANGLE_H
#define __rcRECTANGLE_H

#include "pair.hpp"
#include <deque>


using namespace std;


// Rectangle Class
// There are 2 uncoverntional operator overloads: & and | for intersect and enclose member functions //


template <class T>
class tRectangle
{
public:
    typedef T value_t;
    typedef tpair<T> point_t;
    typedef tRectangle<value_t> self_t;
    

    tRectangle();
    /*
         Makes this rectangle have origin (0,0), with width and
	       m_height both equal to zero.
  */
    tRectangle(T w, T h);
    /*
         Makes this rectangle have the indicated width and
	       m_height.  If w and h are both non-negative, the origin
	       will be (0,0).
    note       If w is negative, the rectangle's origin x-component
	       will be w (negative), and its width will be -w
	       (positive).
	       Similarly, if h is negative, the rectangle's origin
	       y-component will be h (negative), and its m_height will
	       be -h (positive).
  */
    tRectangle(T x, T y, T w, T h);
    /*
         Makes this rectangle have the indicated origin, width,
	       and m_height.
    note       If w is negative, the rectangle's origin x-component
	       will be x+w, and its width will be -w (positive).
	       Similarly, if h is negative, the rectangle's origin
	       y-component will be y+h, and its m_height will be -h
	       (positive).
  */
    tRectangle(const tpair<T> & v1, const tpair<T> & v2);
    /*
         Makes this rectangle the minimum enclosing rectangle
	       of the indicated two points.
  */

    template<typename U>
    tRectangle<U> cast_value_type_to() const;
    /*
        Return a copy of this rectangle in the given type
     */
    

    /* default copy ctor, assignment, dtor OK */

    bool operator==(const tRectangle<T> &) const;
    bool operator!=(const tRectangle<T> & rhs) const { return !(*this == rhs); }
    /*
         Two rectanges are equal iff they have the same origin,
	       width, and m_height.
  */

    const tpair<T> & origin() const { return mUpperLeft; }
    /*
         Returns the origin of this rectangle.  This is its
	       upper-left corner.
  */
    void origin(const tpair<T> &);
    /*
         Sets the origin of this rectangle without modifying
	       its width or m_height.
  */
    void translate(const tpair<T> &);
    /*
         Moves the origin of this rectangle without modifying
	       its width or m_height, by adding the indicated vector
	       displacement to this rectangle's origin.
  */

    const tpair<T> & ul() const { return mUpperLeft; } /* same as origin() */
    tpair<T> ur() const;
    tpair<T> ll() const;
    const tpair<T> & lr() const { return mLowerRight; }
    /*
         Returns the coordinates of the indicated corner of
	       this rectangle.
  */

    tpair<T> center() const { return ((ul() + lr()) / (T)2); }
    /*
         Returns the coordinates of the simple center of
	       this rectangle.
  */
    T width() const;
    /*
         Returns the width of this rectangle.
  */
    void width(T w);
    /*
         Sets the width of this rectangle without changing its
               origin.
    note       If w is negative, its value is added to the x
	       component of the rectangle's origin, and the
	       rectangle's width is set to -w (positive).
  */
    T m_height() const;
    /*
         Returns the m_height of this rectangle.
  */
    void m_height(T h);
    /*
         Sets the m_height of this rectangle without changing its
               origin.
    note       If h is negative, its value is added to the y
	       component of the rectangle's origin, and the
	       rectangle's m_height is set to -h (positive).
  */

    T area() const;
    /*
         Returns the area: width * m_height of this rect
  */

    tpair<T> size() const;
    /*
         Returns the size of this rectangle (the width is the x
	       component of the returned pair; the m_height is the y
	       component).
  */
    void size(const tpair<T> &);
    /*
         Sets the width from the x component of the argument,
	       and sets the m_height from the y component of the
	       argument.
    note       See width(T) and m_height(T) for a
	       description of how negative size values are handled.
  */

    bool isNull(void) const;
    /*
         Tells if this rectangle is a null rectangle (i.e., m_height==0
		OR width==0).
  */

    bool overlaps(const tRectangle<T> &, bool touching = false) const;
    /*
         Tells if the indicated rectangle overlaps this
	       rectangle.
	       If the touching flag is set, rectangles that just
	       barely touch at their borders are considered to
	       overlap.
	       If the flag is not set, rectangles are considered to
	       overlap only if a portion of one rectangle's border is
	       within the other rectangle.
  */

    bool contains(const tRectangle<T> &) const;
    /*
         Tells if the indicated rectangle is contained within
	       this rectangle.
  */

    bool contains(const int, const int) const;
    /*
         Tells if the indicated point is contained within
	       this rectangle.
  */
    void include (const tpair<T> &);
    void include (T x, T y );
    /*
          modify this rectangle to include this point
     
     */
    
    bool contains(const tpair<T> &) const;
    /*
         Tells if the indicated point is contained within
	       this rectangle.
  */

    tRectangle<T> intersect(const tRectangle<T> & other,
                            bool & success) const;
    /*
         Returns a rectangle which is the intersection of the
	       two rectangles.
	       If the two rectangles do not overlap (as indicated
	       by overlaps(other, true)),  and throwIfNoOverlap
	       is false, then a null rectangle is returned; the
	       returned rectangle's origin is undefined; it is the
	       client's responsibility to not use the null rectangle
	       returned by this function for further operations.

    note       Rectangles that just barely touch will produce a null
	       Rectangle whose origin is valid.
  */
    tRectangle<T> operator&(const tRectangle<T> & other) const;
    /*
         tRectangle<T>::intersect(other, false).
         Greatest Common Rectangle of several rectangles can be
	       computed as destRect = (rect1 & rect2 & rect3 & rect4 & ...);
	 followed by a check for a null rectangle
  */
    tRectangle<T> & operator&=(const tRectangle<T> & other);
    /*
         Modifies this rectangle to be the intersection (operator&)
               of the two rectangles. Returns a reference to this object.
  */

    tRectangle<T> enclose(const tRectangle<T> &) const;

    tRectangle<T> operator|(const tRectangle<T> &) const;
    /*
         Returns a rectangle which is the minimum enclosing
	       rectangle of the two rectangles.
	       If tRectangle<T>s A or B are null, (A|B) will be the same as
	       other (non-null) tRectangle<T>.
	       If tRectangle<T>s A and B are both null, (A|B) will return a
	       null tRectangle<T> whose origin is undefined, unless,
	       A.origin()==B.origin().
    note
	       The minimum enclosing rectangle of several rectangles can be
	       computed compactly as
	          destRect = (rect1 | rect2 | rect3 | rect4 | ...);
  */
    tRectangle<T> & operator|=(const tRectangle<T> &);
    /*
         Modifies this rectangle to be the minimum enclosing
	       rectangle of the two rectangles.
  */

    tRectangle<T> trim(T left, T right, T top, T bottom, bool & success) const;

    tRectangle<T> trim(T pad, bool & success) const;

    /*
         Returns a rectangle which is a border adjusted copy of this
               rectangle. The arguments indicate how much the returned
               rectangle is to be shrunk at each of its borders.
               Negative argument values are permitted and mean
               that the rectangle is to be grown rather than shrunk
               at the indicated border.
               If left+right > width() or top+bottom > m_height(),
               and throwIfNegativeSize is false, then negative size
               values are handled as described in width(T)
               and m_height(T).
  */

    tRectangle<T> transpose() const;
    /*
         Returns the transpose of this object. The transpose
	       operation is defined as follows:
	       If two tRectangle<T>s A and B are related by A = B.transpose(), then
		   	A.origin().x() = B.origin.y();
		   	A.origin().y() = B.origin.x();
		   	A.m_height()     = B.width();
		   	A.width()      = B.m_height();
  */

    friend ostream & operator<<(ostream & ous, const tRectangle<T> & dis)
    {
        ous << dis.ul() << "," << dis.lr();
        return ous;
    }


protected:
    void normalize(); /* fix negative width and/or m_height */

    tpair<T> mUpperLeft;
    tpair<T> mLowerRight;
};


// .tamplate keyword make cast_value_type_to a dependent tempalate type
template<class T>
template<typename U>
tRectangle<U> tRectangle<T>::cast_value_type_to() const
{
   return tRectangle<U> ( mUpperLeft.template cast_value_type_to<U>(), mLowerRight.template cast_value_type_to<U>() );
}

template<class T>
inline void tRectangle<T>::include(const tpair<T>& point)
{
    include(point.x(), point.y());
}

template<class T>
inline void tRectangle<T>::include(T x, T y)
{
    if (mUpperLeft.x() > x) mUpperLeft.first = x;
    if (mLowerRight.x() < x) mLowerRight.first = x;
    if (mUpperLeft.y() > y) mUpperLeft.second = y;
    if (mLowerRight.y() < y) mLowerRight.second = y;
    
}

template <class T>
inline tpair<T> tRectangle<T>::ur() const
{
    return tpair<T>(mLowerRight.x(), mUpperLeft.y());
}

template <class T>
inline tpair<T> tRectangle<T>::ll() const
{
    return tpair<T>(mUpperLeft.x(), mLowerRight.y());
}

template <class T>
inline T tRectangle<T>::width() const
{
    return mLowerRight.x() - mUpperLeft.x();
}

template <class T>
inline T tRectangle<T>::m_height() const
{
    return mLowerRight.y() - mUpperLeft.y();
}

template <class T>
inline tpair<T> tRectangle<T>::size() const
{
    return tpair<T>(width(), m_height());
}

template <class T>
inline T tRectangle<T>::area() const
{
    return width() * m_height();
}

template <class T>
inline bool tRectangle<T>::isNull(void) const
{
    return (width() == 0 || m_height() == 0);
}

template <class T>
inline tRectangle<T> tRectangle<T>::operator&(const tRectangle<T> & other) const
{
    return intersect(other, false);
}

template <class T>
inline tRectangle<T> tRectangle<T>::operator|(const tRectangle<T> & other) const
{
    return enclose(other);
}

typedef tRectangle<float> fRect;
typedef tRectangle<double> dRect;
typedef tRectangle<int32_t> iRect;
typedef tRectangle<uint32_t> uiRect;


template <class T>
tRectangle<T>::tRectangle()
    : mUpperLeft(0, 0), mLowerRight(0, 0)
{
}


// We always call normalize. Normalize for unsigned type does nothing. This is a bit inefficient.
template <class T>
tRectangle<T>::tRectangle(T w, T h)
    : mUpperLeft(0, 0), mLowerRight(w, h)
{

    normalize();
}

template <class T>
tRectangle<T>::tRectangle(T x, T y, T w, T h)
    : mUpperLeft(x, y), mLowerRight(x + w, y + h)
{
    normalize();
}

template <class T>
tRectangle<T>::tRectangle(const tpair<T> & v1, const tpair<T> & v2)
    : mUpperLeft(std::min(v1.x(), v2.x()), std::min(v1.y(), v2.y())), mLowerRight(std::max(v1.x(), v2.x()), std::max(v1.y(), v2.y()))
{
    normalize ();
}

template <class T>
void tRectangle<T>::normalize()
{
    T w = width();
    T h = m_height();

    if (w < 0 || h < 0)
    {
        T x = mUpperLeft.x();
        T y = mUpperLeft.y();
        if (w < 0)
        {
            x += w;
            w = -w;
        }
        if (h < 0)
        {
            y += h;
            h = -h;
        }

        mUpperLeft = tpair<T>(x, y);
        mLowerRight = tpair<T>(x + w, y + h);
    }
}

template <class T>
bool tRectangle<T>::operator==(const tRectangle<T> & rhs) const
{
    return (mUpperLeft.x() == rhs.mUpperLeft.x()) &&
           (mUpperLeft.y() == rhs.mUpperLeft.y()) &&
           (mLowerRight.x() == rhs.mLowerRight.x()) &&
           (mLowerRight.y() == rhs.mLowerRight.y());
}

template <class T>
void tRectangle<T>::origin(const tpair<T> & newOrigin)
{
    tpair<T> temp(newOrigin.x() - mUpperLeft.x(),
                  newOrigin.y() - mUpperLeft.y());
    mUpperLeft = newOrigin;
    mLowerRight.x() += temp.x();
    mLowerRight.y() += temp.y();
}
template <class T>
void tRectangle<T>::translate(const tpair<T> & disp)
{
    mUpperLeft.x() += disp.x();
    mUpperLeft.y() += disp.y();
    mLowerRight.x() += disp.x();
    mLowerRight.y() += disp.y();
}

template <class T>
void tRectangle<T>::width(T newWidth)
{
    mLowerRight.x() = mUpperLeft.x() + newWidth;
    normalize();
}

template <class T>
void tRectangle<T>::m_height(T newm_height)
{
    mLowerRight.y() = mUpperLeft.y() + newm_height;
    normalize();
}

template <class T>
void tRectangle<T>::size(const tpair<T> & newSize)
{
    mLowerRight.x() = mUpperLeft.x() + newSize.x();
    mLowerRight.y() = mUpperLeft.y() + newSize.y();
    normalize();
}

template <class T>
bool tRectangle<T>::overlaps(const tRectangle<T> & other, bool touching) const
{
    tpair<T> overlapUpperLeft(std::max(mUpperLeft.x(), other.mUpperLeft.x()),
                              std::max(mUpperLeft.y(), other.mUpperLeft.y()));
    tpair<T> overlapLowerRight(std::min(mLowerRight.x(), other.mLowerRight.x()),
                               std::min(mLowerRight.y(), other.mLowerRight.y()));

    if (touching)
    {
        return (overlapUpperLeft.x() <= overlapLowerRight.x()) &&
               (overlapUpperLeft.y() <= overlapLowerRight.y());
    }
    else
    {
        /* must actually overlap */
        return (overlapUpperLeft.x() < overlapLowerRight.x()) &&
               (overlapUpperLeft.y() < overlapLowerRight.y());
    }
}

template <class T>
bool tRectangle<T>::contains(const tRectangle<T> & other) const
{
    return (other.ul().x() >= ul().x()) &&
           (other.lr().x() <= lr().x()) &&
           (other.ul().y() >= ul().y()) &&
           (other.lr().y() <= lr().y());
}

template <class T>
bool tRectangle<T>::contains(const tpair<T> & point) const
{
    return ul().x() <= point.x && point.x() < lr().x() &&
           ul().y() <= point.y() && point.y() < lr().y();
}


template <class T>
bool tRectangle<T>::contains(const int point_x, const int point_y) const
{
    return ul().x() <= point_x && point_x < lr().x() &&
           ul().y() <= point_y && point_y < lr().y();
}


template <class T>
tRectangle<T> tRectangle<T>::intersect(const tRectangle<T> & other,
                                       bool & success) const
{
    tpair<T> overlapUpperLeft(std::max(mUpperLeft.x(), other.mUpperLeft.x()),
                              std::max(mUpperLeft.y(), other.mUpperLeft.y()));
    tpair<T> overlapLowerRight(std::min(mLowerRight.x(), other.mLowerRight.x()),
                               std::min(mLowerRight.y(), other.mLowerRight.y()));

    success = overlapUpperLeft.x() <= overlapLowerRight.x() && overlapUpperLeft.y() <= overlapLowerRight.y();
    if (!success)
        overlapLowerRight = overlapUpperLeft; // Force a null rectangle

    return tRectangle<T>(overlapUpperLeft, overlapLowerRight);
}

template <class T>
tRectangle<T> & tRectangle<T>::operator&=(const tRectangle<T> & other)
{
    mUpperLeft.x() = std::max(mUpperLeft.x(), other.mUpperLeft.x());
    mUpperLeft.y() = std::max(mUpperLeft.y(), other.mUpperLeft.y());
    mLowerRight.x() = std::min(mLowerRight.x(), other.mLowerRight.x());
    mLowerRight.y() = std::min(mLowerRight.y(), other.mLowerRight.y());

    if ((mUpperLeft.x() > mLowerRight.x()) ||
        (mUpperLeft.y() > mLowerRight.y()))
        mLowerRight = mUpperLeft; // Force a null rectangle
    return *this;
}

template <class T>
tRectangle<T> tRectangle<T>::enclose(const tRectangle<T> & other) const
{
    if (other.isNull())
        return *this;
    if (isNull())
        return other;
    tpair<T> tempMin(std::min(mUpperLeft.x(), other.mUpperLeft.x()),
                     std::min(mUpperLeft.y(), other.mUpperLeft.y()));
    tpair<T> tempMax(std::max(mLowerRight.x(), other.mLowerRight.x()),
                     std::max(mLowerRight.y(), other.mLowerRight.y()));

    return tRectangle<T>(tempMin, tempMax);
}

template <class T>
tRectangle<T> & tRectangle<T>::operator|=(const tRectangle<T> & other)
{
    if (!other.isNull())
    {
        if (isNull())
            *this = other;
        else
        {
            mUpperLeft.x() = std::min(mUpperLeft.x(), other.mUpperLeft.x());
            mUpperLeft.y() = std::min(mUpperLeft.y(), other.mUpperLeft.y());
            mLowerRight.x() = std::max(mLowerRight.x(), other.mLowerRight.x());
            mLowerRight.y() = std::max(mLowerRight.y(), other.mLowerRight.y());
        }
    }

    return *this;
}


template <class T>
tRectangle<T> tRectangle<T>::trim(T pad, bool & success) const
{
    return trim(pad, pad, pad, pad, success);
}

template <class T>
tRectangle<T> tRectangle<T>::trim(T left, T right, T top,
                                  T bottom, bool & success) const
{
    T newWidth = width() - (left + right);
    T newm_height = m_height() - (top + bottom);
    success = newWidth > 0 && newm_height > 0;
    if (success)
        return tRectangle<T>(mUpperLeft + tpair<T>(left, top),
                             mLowerRight - tpair<T>(right, bottom));
    return *this;
}

template <class T>
tRectangle<T> tRectangle<T>::transpose() const
{
    tpair<T> transposeUpperLeft(mUpperLeft.y(), mUpperLeft.x());
    tpair<T> transposeLowerRight(mLowerRight.y(), mLowerRight.x());
    return tRectangle<T>(transposeUpperLeft, transposeLowerRight);
}


template <>
void tRectangle<uint32_t>::normalize();
template <>
tRectangle<uint32_t> tRectangle<uint32_t>::trim(uint32_t left, uint32_t right, uint32_t top, uint32_t bottom, bool &) const;



#endif
