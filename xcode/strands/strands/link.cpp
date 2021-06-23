
#include <vector>
#include <algorithm>
#include "link.h"
#include "contour.h"
#include "chord.h"
#include "region.h"
#include "utils.hpp"
#include "strands.h"
#include <cmath>
#include <limits>
#include <numeric>
#include <math.h>

using namespace std;
using namespace stegers;




std::vector<std::vector<std::vector<int>>> dirtab  {
     { {  1, 0 }, {  1,-1 }, {  1, 1 } },
     { {  1, 1 }, {  1, 0 }, {  0, 1 } },
     { {  0, 1 }, {  1, 1 }, { -1, 1 } },
     { { -1, 1 }, {  0, 1 }, { -1, 0 } },
     { { -1, 0 }, { -1, 1 }, { -1,-1 } },
     { { -1,-1 }, { -1, 0 }, {  0,-1 } },
     { {  0,-1 }, { -1,-1 }, {  1,-1 } },
     { {  1,-1 }, {  0,-1 }, {  1, 0 } }
};

std::vector<std::vector<std::vector<int>>> cleartab = {
  { {  0, 1 }, {  0,-1 } },
  { { -1, 1 }, {  1,-1 } },
  { { -1, 0 }, {  1, 0 } },
  { { -1,-1 }, {  1, 1 } },
  { {  0,-1 }, {  0, 1 } },
  { {  1,-1 }, { -1, 1 } },
  { {  1, 0 }, { -1, 0 } },
  { {  1, 1 }, { -1,-1 } }
};



// eigenvalvect[0][0] - first eigenvector
// eigenvalvect[0][1] - second eigenvector
// eigenvalvect[1][0] - first eigenvector x
// eigenvalvect[1][1] - first eigenvector y
// eigenvalvect[2][0] - second eigenvector x
// eigenvalvect[2][1] - second eigenvector y

// Compute the eigenvalues and eigenvectors of the Hessian matrix given by dfdrr, dfdrc, and dfdcc,
// and sort them in descending order according to their absolute values

void compute_eigenvals(double dfdrr, double dfdrc,
    double dfdcc, std::vector<std::vector<double>>& eigenvalvect) {

    double theta, t, n1, n2; // , phi;
    double c = 1.0;
    double s = 0.0;
    double e1 = dfdrr;
    double e2 = dfdcc;

    /* Compute the eigenvalues and eigenvectors of the Hessian matrix. */
    if (dfdrc != 0.0) {
        theta = 0.5 * (dfdcc - dfdrr) / dfdrc;
        t = 1.0 / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        if (theta < 0.0) {
            t = -t;
        }
        c = 1.0 / std::sqrt(t * t + 1.0);
        s = t * c;
        e1 = dfdrr - t * dfdrc;
        e2 = dfdcc + t * dfdrc;
    }

    n1 = c;
    n2 = -s;

    /* If the absolute value of an eigenvalue is larger than the other, put that
       eigenvalue into first position.  If both are of equal absolute value, put
       the negative one first. */
    if (std::abs(e1) > std::abs(e2)) {
        eigenvalvect[0][0] = e1;
        eigenvalvect[0][1] = e2;
        eigenvalvect[1][0] = n1;
        eigenvalvect[1][1] = n2;
        eigenvalvect[2][0] = -n2;
        eigenvalvect[2][1] = n1;
    }
    else if (std::abs(e1) < std::abs(e2)) {
        eigenvalvect[0][0] = e2;
        eigenvalvect[0][1] = e1;
        eigenvalvect[1][0] = -n2;
        eigenvalvect[1][1] = n1;
        eigenvalvect[2][0] = n1;
        eigenvalvect[2][1] = n2;
    }
    else {
        if (e1 < e2) {
            eigenvalvect[0][0] = e1;
            eigenvalvect[0][1] = e2;
            eigenvalvect[1][0] = n1;
            eigenvalvect[1][1] = n2;
            eigenvalvect[2][0] = -n2;
            eigenvalvect[2][1] = n1;
        }
        else {
            eigenvalvect[0][0] = e2;
            eigenvalvect[0][1] = e1;
            eigenvalvect[1][0] = -n2;
            eigenvalvect[1][1] = n1;
            eigenvalvect[2][0] = n1;
            eigenvalvect[2][1] = n2;
        }
    }

}


// store results of eigendecomposition
    //	[0] = first eigenvalue
    //	[1] = eigenvector coordinate y
    //	[2] = eigenvector coordinate x
    //	[3] = derivative of image in y
    //	[4] = derivative of image in x
    //	[5] = "super-resolved" y		t_y, or dlpy
    //	[6] = "super-resolved" x		t_x, or dlpx	

void link::compute_line_points(const std::vector<std::vector<float>>& derivs,
    std::vector<int32_t>& ismax, std::vector<std::vector<double>>& ennpp,
    double low, double high) {

    int r, c;
    std::vector<std::vector<double>> line_points = RectangularDoubleVector(7, static_cast<int>(total()));
    std::vector<std::vector<double>> eigenvalvect = RectangularDoubleVector(3, 2);
    ismax.resize(total(), 0);

    timer::precise_stopwatch stopwatch;

    //all derivatives
    int l;
    double a, b, n1, n2, t;

    /*compute eigenval and vectors */
    for (r = 0; r < m_height; r++) {
        float* row_ptr = (float*)(m_espace.ptr(r));
        for (c = 0; c < m_width; c++, row_ptr++) {
            l = r * m_width + c;
            compute_eigenvals(derivs[2][l], derivs[3][l], derivs[4][l], eigenvalvect);
            auto val = m_polarity == strands::line_polarity::light ? -eigenvalvect[0][0] : eigenvalvect[0][0];

            if (val > 0.0) {
                *row_ptr = (float)val;
                line_points[0][l] = val; // first eigenvalue
                n1 = eigenvalvect[1][0]; // igenvector coordinate y
                n2 = eigenvalvect[1][1]; // eigenvector coordinate x
                a = derivs[2][l] * n1 * n1 + 2.0 * derivs[3][l] * n1 * n2 + derivs[4][l] * n2 * n2;
                b = derivs[0][l] * n1 + derivs[1][l] * n2;
                if (a == 0.0) continue;
                t = (-1) * b / a;
                auto p1 = t * n1;
                auto p2 = t * n2;
                if (std::fabs(p1) > 0.6 || std::fabs(p2) > 0.6) continue;
                ismax[l] = val >= high ? 2 : val >= low ? 1 : 0;
                line_points[5][l] = r + p1;
                line_points[6][l] = c + p2;
                line_points[1][l] = n1; //  "super-resolved" y
                line_points[2][l] = n2; //  "super-resolved" x
                line_points[3][l] = derivs[0][l];
                line_points[4][l] = derivs[1][l];

            }
        }
    }
    ennpp = line_points;

    if (debug()) {
        auto actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        std::cout << "Computing Line Points Time " << actual_wait_time << " microseconds" << std::endl;;
    }
}

void link::getismax(double low, double high, const vector<double>& eigval, const vector<double>& posx, const vector<double>& posy, std::vector<uint32_t>& iismax){

    static int32_t zero = 0;
    iismax.resize(total (), zero);
    int l;
    double val;
    for (int py = 0; py < height(); ++py) {
        for (int px = 0; px < width(); ++px) {
            l = LCOR(py, px, width());
            val = eigval[l];
            if (!isnan(posx[l]) && std::fabs(posy[l] - px) <= m_pixel_boundary && std::fabs(posx[l] - py) <= m_pixel_boundary && val >= low)
                iismax[l] = (val >= high) ? 2 : 1;
        }
    }
}

void link::interpolate_gradient(const std::vector<double>& gradx, const std::vector<double>& grady, double px, double py, doublepoint& gradxy)
{
    long   gix, giy, gpos;
    double gfx, gfy, gx1, gy1, gx2, gy2, gx3, gy3, gx4, gy4;

    gix = (long)std::floor(px);
    giy = (long)std::floor(py);
    gfx = px; // % 1.0; //check whether it works as promised
    gfy = py; // % 1.0;
    gpos = LCOR(gix, giy, width());
    gx1 = gradx[(int)gpos];
    gy1 = grady[(int)gpos];
    gpos = LCOR(gix + 1, giy, width());
    gx2 = gradx[(int)gpos];
    gy2 = grady[(int)gpos];
    gpos = LCOR(gix, giy + 1, width());
    gx3 = gradx[(int)gpos];
    gy3 = grady[(int)gpos];
    gpos = LCOR(gix + 1, giy + 1, width());
    gx4 = gradx[(int)gpos];
    gy4 = grady[(int)gpos];
    gradxy.cx = (1 - gfy) * ((1 - gfx) * gx1 + gfx * gx2) + gfy * ((1 - gfx) * gx3 + gfx * gx4);
    gradxy.cy = (1 - gfy) * ((1 - gfx) * gy1 + gfx * gy2) + gfy * ((1 - gfx) * gy3 + gfx * gy4);
}

double link::interpolate_response(const std::vector<double>& resp, int32_t x, int32_t y, double px, double py)
{


    double i1, i2, i3, i4, i5, i6, i7, i8, i9;
    double t1, t2, t3, t4, t5, t6;
    double d, dr, dc, drr, drc, dcc;
    double xx, yy;

    i1 = resp[LCOR(BR(x - 1), BC(y - 1),width())];
    i2 = resp[LCOR(BR(x - 1), y,width())];
    i3 = resp[LCOR(BR(x - 1), BC(y + 1),width())];
    i4 = resp[LCOR(x, BC(y - 1),width())];
    i5 = resp[LCOR(x, y,width())];
    i6 = resp[LCOR(x, BC(y + 1),width())];
    i7 = resp[LCOR(BR(x + 1), BC(y - 1),width())];
    i8 = resp[LCOR(BR(x + 1), y,width())];
    i9 = resp[LCOR(BR(x + 1), BC(y + 1),width())];
    t1 = i1 + i2 + i3;
    t2 = i4 + i5 + i6;
    t3 = i7 + i8 + i9;
    t4 = i1 + i4 + i7;
    t5 = i2 + i5 + i8;
    t6 = i3 + i6 + i9;
    d = (-i1 + 2 * i2 - i3 + 2 * i4 + 5 * i5 + 2 * i6 - i7 + 2 * i8 - i9) / 9;
    dr = (t3 - t1) / 6;
    dc = (t6 - t4) / 6;
    drr = (t1 - 2 * t2 + t3) / 6;
    dcc = (t4 - 2 * t5 + t6) / 6;
    drc = (i1 - i3 - i7 + i9) / 4;
    xx = px - x;
    yy = py - y;
    return d + xx * dr + yy * dc + xx * xx * drr + xx * yy * drc + yy * yy * dcc;
}


void link::closest_point(double lx, double ly, double dx, double dy, double px, double py, doublepoint& out){

    double mx, my, den, nom, tt;

    mx = px - lx;
    my = py - ly;
    den = dx * dx + dy * dy;
    nom = mx * dx + my * dy;
    if (den != 0)
        tt = nom / den;
    else
        tt = 0;
    out.cx = lx + tt * dx;
    out.cy = ly + tt * dy;
    out.t = tt;
}


void link::bresenham(double nx, double ny, double px, double py, double length, std::vector<offset>& out) {
    {

        out.resize(0);
        int i, x, y, s1, s2, xchg, maxit;
        double e, dx, dy, t;

        x = 0;
        y = 0;
        dx = std::fabs(nx);
		dy = std::fabs(ny);
        s1 = sgn(nx);
        s2 = sgn(ny);
        px *= s1;
        py *= s2;
        if (dy > dx) {
            t = dx;
            dx = dy;
            dy = t;
            t = px;
            px = py;
            py = t;
            xchg = 1;
        }
        else {
            xchg = 0;
        }
        maxit = (int)std::ceil(length * dx);
        e = (0.5 - px) * dy / dx - (0.5 - py);
        for (i = 0; i <= maxit; i++) {
            out.emplace_back(x, y);

            while (e >= -1e-8) {
                if (std::fabs(xchg) > 0) x += s1;
                else y += s2;
                e--;
                if (e > -1) {
                    out.emplace_back(x, y);

                }
            }
            if (std::fabs(xchg) > 0) y += s2;
            else x += s1;
            e += dy / dx;
        }
    }
}


void link::compute_contours(std::vector<int32_t> ismax, const std::vector<double>& eigval,
    const std::vector<double>& posx, const std::vector<double>& posy,
    const std::vector<double>& normx, const std::vector<double>& normy,
    const std::vector<double>& gradx, const std::vector<double>& grady){

    long area;
    int i, k, l, pos, nexti;
    int nextpos;
    int it;
    long num_pnt, num_cont, num_junc;
    long x, y;
    long begin, end;
    long indx_max;
    double max;
    long maxx, maxy, nextx, nexty;
    double nx, ny, mx, my;
    double px, py, nextpx, nextpy;
    double alpha, nextalpha, beta, last_beta;
    int octant, last_octant;
    bool nextismax;
    double diff, mindiff, diff1, diff2, dist, mindist;
    double dx, dy;
    double s, t, gx, gy;
    double length, response;
    int num_add;
    int m = 0;
    int j = 0;
    double end_angle = 0;
    double end_resp = 0;
    contour tmp_cont;
    bool add_ext;
    doublepoint closestpnt;
    double MAX_ANGLE_DIFFERENCE = pi / 6.0;
    std::vector<contour> cont;
    timer::precise_stopwatch stopwatch;


    contour_class cls;
    vector<float> row ;
    vector<float> col ;
    vector<float> angle ;
    vector<float> resp ;


    vector<float> extx;
    vector<float> exty;
    vector<offset> line;
    vector<float> trow;
    vector<float> tcol;
    vector<float> tangle;
    vector<float> tresp;

    vector<chord> rl;
    vector<short> label(total(), 0);
    vector<int32_t> indx(total(), 0);


    // Select all pixels that can be starting points for lines.           
    region seg (ismax, 2, width(), height());

    // Count the number of possible starting points. 
    area = 0;
    for (i = 0; i < seg.num; i++) {
        area += seg.rl[i].ce - seg.rl[i].cb + 1;
    }

    // Create the index of possible starting points. 
    std::vector<crossRef> cross;

    rl = seg.rl;
    for (i = 0; i < seg.num; i++)
    {
        x = rl[i].r;
        for (y = rl[i].cb; y <= rl[i].ce; y++)
        {
            pos = (int)LCOR(x, y, width());
            cross.emplace_back((short)x, (short)y, std::fabs(eigval[pos]), false);
        }
    }

    sort(cross.begin(), cross.end(), greater<crossRef> ());

    if (area > m_max_line_pts) {
        m_crowded = true;
        area = m_max_line_pts;
    }

     for (i = 0; i < area; i++)
        indx[LCOR(cross[i].x(), cross[i].y(), width())] = i + 1;

    //reset everything
    num_cont = 0;
    num_junc = 0;
    vector<junction>junc;

    // Link lines points. 
    indx_max = 0;
    for (;;)
    {
        // Contour class unknown at this point; therefore assume both ends free. 
        cls = contour_class::cont_no_junc;
        while (indx_max < area && cross[indx_max].done())
            indx_max++;
        // Stop if no feasible starting point exists. 
        if (indx_max == area)
            break;
        max = cross[indx_max].value();
        maxx = cross[indx_max].x();
        maxy = cross[indx_max].y();
        if (max == 0.0)
            break;
        row.resize(0) ;
        col.resize(0);
        resp.resize(0);
        angle.resize(0);

        // Add starting point to the line. 
        num_pnt = 0;
        pos = (int)LCOR(maxx, maxy, width());
        label[pos] = (short)(num_cont + 1);
        if ( indx[pos] != 0)
            cross[indx[pos] - 1].setDone();
        row.push_back((float)posx[pos]);
        col.push_back((float)posy[pos]);

        // Select line direction. 
        nx = -normy[pos];
        ny = normx[pos];
        alpha = std::atan2(ny, nx);
        if (alpha < 0.0)
            alpha += 2.0 * pi;
        if (alpha >= pi)
            alpha -= pi;
        octant = (int)(std::floor(4.0 / pi * alpha + 0.5)) % 4;
        // Select normal to the line.  The normal points to the right of the line
        // as the line is traversed from 0 to num-1.  Since the points are sorted
        // in reverse order before the second iteration, the first beta actually
        // has to point to the left of the line! 
        beta = alpha + pi / 2.0;
        if (beta >= 2.0 * pi)
            beta -= 2.0 * pi;
        angle.push_back((float)beta);
        resp.push_back((float)interpolate_response(eigval, maxx, maxy, posx[pos], posy[pos]));
        num_pnt++;

        // Mark double responses as processed. 
        for (i = 0; i < 2; i++)
        {
            nextx = maxx + cleartab[octant][i][0];
            nexty = maxy + cleartab[octant][i][1];
            if (nextx < 0 || nextx >= height() || nexty < 0 || nexty >= width())
                continue;
            nextpos = (int)LCOR(nextx, nexty, width());
            if (ismax[nextpos] > 0)
            {
                nextalpha = get_angle(normx, normy, nextpos);
                diff = std::abs(alpha - nextalpha);
                if (diff >= pi / 2.0)
                    diff = pi - diff;
                if (diff < MAX_ANGLE_DIFFERENCE)
                {
                    label[nextpos] = (short)(num_cont + 1);
                    if ( indx[nextpos] != 0)
                        cross[indx[nextpos] - 1].setDone();
                }
            }
        }

        for (it = 1; it <= 2; it++)
        {
            if (it == 1)
            {
                // Search along the initial line direction in the first iteration. 
                x = maxx;
                y = maxy;
                pos = (int)LCOR(x, y, width());
                alpha = get_angle(normx, normy, pos);
                last_octant = (int)(std::floor(4.0 / pi * alpha + 0.5)) % 4;
                last_beta = alpha + pi / 2.0;
                if (last_beta >= 2.0 * pi)
                    last_beta -= 2.0 * pi;
            }
            else
            {
                // Search in the opposite direction in the second iteration. 
                x = maxx;
                y = maxy;
                pos = (int)LCOR(x, y, width());
                alpha = get_angle(normx, normy, pos);
                last_octant = (int)(std::floor(4.0 / pi * alpha + 0.5)) % 4 + 4;
                last_beta = alpha + pi / 2.0;
                if (last_beta >= 2.0 * pi)
                    last_beta -= 2.0 * pi;
            }
            if (it == 2)
            {
                // Sort the points found in the first iteration in reverse. 
                std::reverse(row.begin(), row.end());
                std::reverse(col.begin(), col.end());
                std::reverse(angle.begin(), angle.end());
                std::reverse(resp.begin(), resp.end());
            }


            // Now start adding appropriate neighbors to the line. 
            for (;;) {
                pos = (int)LCOR(x, y, width());
                px = posx[pos];
                py = posy[pos];
                // Orient line direction w.r.t. the last line direction. 
                alpha = get_angle(normx, normy, pos);
                octant = (int)(std::floor(4.0 / pi * alpha + 0.5)) % 4;
                switch (octant) {
                case 0:
                    if (last_octant >= 3 && last_octant <= 5)
                        octant = 4;
                    break;
                case 1:
                    if (last_octant >= 4 && last_octant <= 6)
                        octant = 5;
                    break;
                case 2:
                    if (last_octant >= 4 && last_octant <= 7)
                        octant = 6;
                    break;
                case 3:
                    if (last_octant == 0 || last_octant >= 6)
                        octant = 7;
                    break;
                }
                last_octant = octant;

                // Determine appropriate neighbor. 
                nextismax = false;
                nexti = 1;
                mindiff = std::numeric_limits<double>::max();
                for (i = 0; i < 3; i++) {
                    nextx = x + dirtab[octant][i][0];
                    nexty = y + dirtab[octant][i][1];
                    if (nextx < 0 || nextx >= height() || nexty < 0 || nexty >= width())
                        continue;
                    nextpos = (int)LCOR(nextx, nexty, width());
                    if (ismax[nextpos] == 0)
                        continue;
                    nextpx = posx[nextpos];
                    nextpy = posy[nextpos];
                    dx = nextpx - px;
                    dy = nextpy - py;
                    dist = std::sqrt(dx * dx + dy * dy);
                    nextalpha = get_angle(normx, normy, nextpos);
                    diff = std::abs(alpha - nextalpha);
                    if (diff >= pi / 2.0)
                        diff = pi - diff;
                    diff = dist + diff;
                    if (diff < mindiff) {
                        mindiff = diff;
                        nexti = i;
                    }
                    if (!ismax[nextpos] == 0)  
						nextismax = true;
                }
                // Mark double responses as processed. 
                for (i = 0; i < 2; i++) {
                    nextx = x + cleartab[octant][i][0];
                    nexty = y + cleartab[octant][i][1];
                    if (nextx < 0 || nextx >= height() || nexty < 0 || nexty >= width())
                        continue;
                    nextpos = (int)LCOR(nextx, nexty, width());
                    if (ismax[nextpos] > 0) {
                        nextalpha = get_angle(normx, normy, nextpos);
                        diff = std::abs(alpha - nextalpha);
                        if (diff >= pi / 2.0)
                            diff = pi - diff;
                        if (diff < MAX_ANGLE_DIFFERENCE) {
                            label[nextpos] = (short)(num_cont + 1);
                            if ( indx[nextpos] !=  0)
                                cross[indx[nextpos] - 1].setDone();
                        }
                    }
                }

                // Have we found the end of the line? 
                if (!nextismax)
                    break;
                // If not, add the neighbor to the line. 
                x += dirtab[octant][nexti][0];
                y += dirtab[octant][nexti][1];

                pos = (int)LCOR(x, y, width());
                row.push_back((float)posx[pos]);
                col.push_back((float)posy[pos]);

                // Orient normal to the line direction w.r.t. the last normal. 
                beta = get_angle(normx, normy, pos);
                diff1 = std::abs(beta - last_beta);
                if (diff1 >= pi)
                    diff1 = 2.0 * pi - diff1;
                diff2 = std::abs(beta + pi - last_beta);
                if (diff2 >= pi)
                    diff2 = 2.0 * pi - diff2;
                if (diff1 < diff2) {
                    angle.push_back((float)beta);
                    last_beta = beta;
                }
                else {
                    angle.push_back((float)(beta + pi));
                    last_beta = beta + pi;
                }

                resp.push_back((float)interpolate_response(eigval, x, y, posx[pos], posy[pos]));
                num_pnt++;

                // If the appropriate neighbor is already processed a junction point is found. 
                if (label[pos] > 0) {

                    // Look for the junction point in the other line. 
                    k = label[pos] - 1;
                    if (k == num_cont) {
                        // Line intersects itself. 
                        for (j = 0; j < num_pnt - 1; j++) {
                            if (row[j] == posx[pos] && col[j] == posy[pos]) {
                                if (j == 0) {
                                    // Contour is closed. 
                                    cls = contour_class::cont_closed;
                                    std::reverse(row.begin(), row.end());
                                    std::reverse(col.begin(), col.end());
                                    std::reverse(angle.begin(), angle.end());
                                    std::reverse(resp.begin(), resp.end());
                                    it = 2;
                                }
                                else {
                                    if (it == 2) {
                                        // Determine contour class. 
                                        if (cls == contour_class::cont_start_junc)
                                            cls = contour_class::cont_both_junc;
                                        else
                                            cls = contour_class::cont_end_junc;
                                        // Index j is the correct index. 
                                        junc.emplace_back(num_cont, num_cont, j, (float)posx[pos], (float)posy[pos]);
                                        num_junc++;
                                    }
                                    else {
                                        // Determine contour class. 
                                        cls = contour_class::cont_start_junc;
                                        // Index num_pnt-1-j is the correct index since the line
                                         //  is going to be sorted in reverse. 
                                        junc.emplace_back(num_cont, num_cont, num_pnt - 1 - j, (float)posx[pos], (float)posy[pos]);
                                        num_junc++;
                                    }
                                }
                                break;
                            }
                        }
                        // Mark this case as being processed for the algorithm below. 
                        j = -1;
                    }
                    else {
                        for (j = 0; j < cont[k].num; j++)
                        {
                            if (cont[k].row[j] == posx[pos] && cont[k].col[j] == posy[pos])
                                break;
                            mindist = std::sqrt(std::pow(cont[k].row[j] - posx[pos], 2) + std::pow(cont[k].col[j] - posy[pos], 2));
                            if (mindist < 0.00001)
                                break;
                        }
                        // If no point can be found on the other line a double response
                        //   must have occured.  In this case, find the nearest point on
                        //   the other line and add it to the current line. 
                        if (j == cont[k].num) {
                            mindist = std::numeric_limits<double>::max();
                            j = -1;
                            for (l = 0; l < cont[k].num; l++) {
                                dx = posx[pos] - cont[k].row[l];
                                dy = posy[pos] - cont[k].col[l];
                                dist = std::sqrt(dx * dx + dy * dy);
                                if (dist < mindist) {
                                    mindist = dist;
                                    j = l;
                                }
                            }
                            // Add the point with index j to the current line. 

                            row.push_back(cont[k].row[j]);
                            col.push_back(cont[k].col[j]);
                            beta = cont[k].angle[j];
                            if (beta >= pi)
                                beta -= pi;
                            diff1 = std::abs(beta - last_beta);
                            if (diff1 >= pi)
                                diff1 = 2.0 * pi - diff1;
                            diff2 = std::abs(beta + pi - last_beta);
                            if (diff2 >= pi)
                                diff2 = 2.0 * pi - diff2;
                            if (diff1 < diff2)
                                angle.push_back((float)beta);
                            else
                                angle.push_back((float)(beta + pi));
                            resp.push_back(cont[k].response[j]);
                            num_pnt++;
                        }
                    }

                    // Add the junction point only if it is not one of the other line's
                    // endpoints. 
                    if (j > 0 && j < cont[k].num - 1) {
                        // Determine contour class. 
                        if (it == 1)
                            cls = contour_class::cont_start_junc;
                        else if (cls == contour_class::cont_start_junc)
                            cls = contour_class::cont_both_junc;
                        else
                            cls = contour_class::cont_end_junc;
                        // Add the new junction. 
                        junc.emplace_back(k, num_cont, j, row[(int)(num_pnt - 1)], col[(int)(num_pnt - 1)]);
                        num_junc++;
                    }
                    break;
                }
                label[pos] = (short)(num_cont + 1);
                if ( indx[pos] !=  0)
                    cross[(int)(indx[pos] - 1)].setDone();
            }
        }


        if (num_pnt > 1) {
            // Only add lines with at least two points.             
            cont.emplace_back(num_pnt, row, col, angle, resp, cls);
            num_cont++;
        }
        else {
            // Delete the point from the label image; we can use maxx and maxy
            // as the coordinates in the label image in this case. 
            for (i = -1; i <= 1; i++) {
                for (j = -1; j <= 1; j++) {
                    pos = (int)LCOR(BR(maxx + i), BC(maxy + j), width());
                    if (label[pos] == num_cont + 1)
                        label[pos] = 0;
                }
            }
        }
    }

    // Now try to extend the lines at their ends to find additional junctions. 
    if (extend_lines()) {
        // Sign by which the gradient has to be multiplied below. 
        s = polarity() == strands::line_polarity::light ? 1 : -1;
        length = MAX_LINE_EXTENSION;

        auto max_line = (int)ceil(length * 3);
        vector<offset> line(max_line);
        std::vector<float> extx(max_line);
        std::vector<float> exty(max_line);


        for (i = 0; i < num_cont; i++) {

            tmp_cont = cont[i];
            num_pnt = tmp_cont.num;
            if (num_pnt == 1)
                continue;
            if (tmp_cont.cont_class == contour_class::cont_closed)
                continue;
            trow = tmp_cont.row;
            tcol = tmp_cont.col;
            tangle = tmp_cont.angle;
            tresp = tmp_cont.response;
            // Check both ends of the line (it==-1: start, it==1: end). 
            for (it = -1; it <= 1; it += 2) {

                // Determine the direction of the search line.  This is done by using 
                  // the normal to the line (angle).  Since this normal may point to
                  // the left of the line (see below) we have to check for this case by
                  // comparing the normal to the direction of the line at its respective
                  // end point. 
                if (it == -1) {
                    // Start point of the line. 
                    if (tmp_cont.cont_class == contour_class::cont_start_junc ||
                        tmp_cont.cont_class == contour_class::cont_both_junc)
                        continue;
                    dx = trow[1] - trow[0];
                    dy = tcol[1] - tcol[0];
                    alpha = tangle[0];
                    nx = std::cos(alpha);
                    ny = std::sin(alpha);
                    if (nx * dy - ny * dx < 0) {
                        // Turn the normal by +90 degrees. 
                        mx = -ny;
                        my = nx;
                    }
                    else {
                        // Turn the normal by -90 degrees. 
                        mx = ny;
                        my = -nx;
                    }
                    px = trow[0];
                    py = tcol[0];
                    response = tresp[0];
                }
                else
                {
                    // End point of the line. 
                    if (tmp_cont.cont_class == contour_class::cont_end_junc ||
                        tmp_cont.cont_class == contour_class::cont_both_junc)
                        continue;
                    dx = trow[(int)(num_pnt - 1)] - trow[(int)(num_pnt - 2)];
                    dy = tcol[(int)(num_pnt - 1)] - tcol[(int)(num_pnt - 2)];
                    alpha = tangle[(int)(num_pnt - 1)];
                    nx = std::cos(alpha);
                    ny = std::sin(alpha);
                    if (nx * dy - ny * dx < 0) {
                        // Turn the normal by -90 degrees. 
                        mx = ny;
                        my = -nx;
                    }
                    else {
                        // Turn the normal by +90 degrees. 
                        mx = -ny;
                        my = nx;
                    }
                    px = trow[(int)(num_pnt - 1)];
                    py = tcol[(int)(num_pnt - 1)];
                    response = tresp[(int)(num_pnt - 1)];
                }
                // Determine the current pixel and calculate the pixels on the search line. 
                x = (long)std::floor(px + 0.5);
                y = (long)std::floor(py + 0.5);
                dx = px - x;
                dy = py - y;
                bresenham(mx, my, dx, dy, length, line);

                // Now determine whether we can go only uphill (bright lines) or
                // downhill (dark lines) until we hit another line. 
                num_add = 0;
                add_ext = false;
                for (k = 0; k < line.size(); k++) {
                    nextx = x + line[k].x;
                    nexty = y + line[k].y;
                    doublepoint closestpnt;
                    closest_point(px, py, mx, my, (double)nextx, (double)nexty, closestpnt);
                    nextpx = closestpnt.cx;
                    nextpy = closestpnt.cy;
                    t = closestpnt.t;

                    // Ignore points before or less than half a pixel away from the
                    // true end point of the line. 
                    if (t <= 0.5)
                        continue;
                    // Stop if the gradient can't be interpolated any more or if the
                    // next point lies outside the image. 
                    if (nextpx < 0 || nextpy < 0 ||
                        nextpx >= height() - 1 || nextpy >= width() - 1 ||
                        nextx < 0 || nexty < 0 ||
                        nextx >= height() || nexty >= width())
                        break;
                    doublepoint interpoint;
                    interpolate_gradient(gradx, grady, nextpx, nextpy,  interpoint);
                    gx = interpoint.cx;
                    gy = interpoint.cy;
                    // Stop if we can't go uphill anymore.  This is determined by the
                    // dot product of the line direction and the gradient.  If it is
                    // smaller than 0 we go downhill (reverse for dark lines). 
                    nextpos = (int)LCOR(nextx, nexty, width());
                    if (s * (mx * gx + my * gy) < 0 && label[nextpos] == 0)
                        break;
                    // Have we hit another line? 
                    if (label[nextpos] > 0) {
                        m = label[nextpos] - 1;
                        // Search for the junction point on the other line. 
                        mindist = std::numeric_limits<double>::max();
                        j = -1;
                        for (l = 0; l < cont[m].num; l++) {
                            dx = nextpx - cont[m].row[l];
                            dy = nextpy - cont[m].col[l];
                            dist = std::sqrt(dx * dx + dy * dy);
                            if (dist < mindist) {
                                mindist = dist;
                                j = l;
                            }
                        }
                        // This should not happen...  But better safe than sorry... 
                        if (mindist > 3.0)
                            break;
                        extx.push_back(cont[m].row[j]);
                        exty.push_back(cont[m].col[j]);
                        end_resp = cont[m].response[j];
                        end_angle = cont[m].angle[j];
                        beta = end_angle;
                        if (beta >= pi)
                            beta -= pi;
                        diff1 = std::abs(beta - alpha);
                        if (diff1 >= pi)
                            diff1 = 2.0 * pi - diff1;
                        diff2 = std::abs(beta + pi - alpha);
                        if (diff2 >= pi)
                            diff2 = 2.0 * pi - diff2;
                        if (diff1 < diff2)
                            end_angle = beta;
                        else
                            end_angle = beta + pi;
                        num_add++;
                   /*     if (DEBUG_show_extensions)
                        {
                            resolved_p = new OvalRoi(cont[m).col[j) - 0.25, cont[m).row[j) - 0.25, 0.5, 0.5);
                            resolved_p.setStrokeColor(Color.GREEN);
                            resolved_p.setPosition(nFrame + 1);
                            image_overlay.add(resolved_p);
                            imp.setOverlay(image_overlay);
                            imp.updateAndRepaintWindow();
                            imp.show();
                        }*/
                        add_ext = true;
                        break;
                    }
                    else {
                        extx.push_back((float)nextpx);
                        exty.push_back((float)nextpy);
                    /*    if (DEBUG_show_extensions)
                        {
                            resolved_p = new OvalRoi(nextpy - 0.25, nextpx - 0.25, 0.5, 0.5);
                            resolved_p.setStrokeColor(Color.GREEN);
                            resolved_p.setPosition(nFrame + 1);
                            image_overlay.add(resolved_p);
                            imp.setOverlay(image_overlay);
                            imp.updateAndRepaintWindow();
                            imp.show();
                        }*/
                        num_add++;
                    }
                }
                if (add_ext) {
                    // Make room for the new points. 
                    num_pnt += num_add;

                    tmp_cont.row = trow;
                    tmp_cont.col = tcol;
                    tmp_cont.angle = tangle;
                    tmp_cont.response = tresp;
                    tmp_cont.num = num_pnt;
                    if (it == -1) {
                        // Move points on the line up num_add places. 

                        // Insert points at the beginning of the line. 
                        for (k = 0; k < num_add; k++)
                        {
                            //cause order of insertion is different
                            trow.insert(trow.begin(), extx[k]);
                            tcol.insert(tcol.begin(), exty[k]);

                            tangle.insert(tangle.begin(), (float)alpha);
                            tresp.insert(tresp.begin(), (float)response);
                        }

                        tangle[0] = (float)end_angle;
                        tresp[0] = (float)end_resp;
                        // Adapt indices of the previously found junctions. 
                        for (k = 0; k < num_junc; k++) {
                            if (junc[k].cont1 == i)
                                junc[k].pos += num_add;
                        }
                    }
                    else {
                        // Insert points at the end of the line. 
                        for (k = 0; k < num_add; k++) {

                            trow.push_back(extx[k]);
                            tcol.push_back(exty[k]);
                            tangle.push_back((float)alpha);
                            tresp.push_back((float)response);
                        }
                        tangle[num_pnt - 1] = (float)end_angle;
                        tresp[num_pnt - 1] = (float)end_resp;
                    }

                    // Add the junction point only if it is not one of the other line's
                     //  endpoints. 
                    if (j > 0 && j < cont[m].num - 1) {
                        if (it == -1) {
                            if (tmp_cont.cont_class == contour_class::cont_end_junc)
                                tmp_cont.cont_class = contour_class::cont_both_junc;
                            else
                                tmp_cont.cont_class = contour_class::cont_start_junc;
                        }
                        else {
                            if (tmp_cont.cont_class == contour_class::cont_start_junc)
                                tmp_cont.cont_class = contour_class::cont_both_junc;
                            else
                                tmp_cont.cont_class = contour_class::cont_end_junc;
                        }
                        junc.emplace_back();
                        junc[(int)num_junc].cont1 = m;
                        junc[(int)num_junc].cont2 = i;
                        junc[(int)num_junc].pos = j;
                        if (it == -1) {
                            junc[(int)num_junc].x = trow[0];
                            junc[(int)num_junc].y = tcol[0];
                        }
                        else {
                            junc[(int)num_junc].x = trow[(int)(num_pnt - 1)];
                            junc[(int)num_junc].y = tcol[(int)(num_pnt - 1)];
                        }
                        num_junc++;
                    }
                }
            }
        }
    }

    // Done with linking.  Now split the lines at the junction points. 

    sort(junc.begin(),junc.end(), greater<junction>());
    m_junctions = junc;

    if (split_lines())
    {

        for (i = 0; i < num_junc; i += k) {
            j = (int)junc[i].cont1;
            tmp_cont = cont[j];
            num_pnt = tmp_cont.num;

            // Count how often line j needs to be split. 
            auto counter = 0;
            for (auto index = 0; index < num_junc; index++) {
                if ((i + k) < num_junc && junc[(i + k)].cont1 == j)
                    counter++;
            }

            if (counter == 1 && tmp_cont.row.size() > (num_pnt - 1) &&
                tmp_cont.row[0] == tmp_cont.row[(int)(num_pnt - 1)] &&
                tmp_cont.col[0] == tmp_cont.col[(int)(num_pnt - 1)])
            {
                // If only one junction point is found and the line is closed it only
                 // needs to be rearranged cyclically, but not split. 
                begin = junc[i].pos;
                trow = tmp_cont.row;
                tcol = tmp_cont.col;
                tangle = tmp_cont.angle;
                tresp = tmp_cont.response;
                for (l = 0; l < num_pnt; l++) {
                    pos = (int)(begin + l);
                    // Skip starting point so that it is not added twice. 
                    if (pos >= num_pnt)
                        pos = (int)(begin + l - num_pnt + 1);
                    tmp_cont.row[l] = trow[pos];
                    tmp_cont.col[l] = tcol[pos];
                    tmp_cont.angle[l] = tangle[pos];
                    tmp_cont.response[l] = tresp[pos];
                }
                // Modify contour class. 
                tmp_cont.cont_class = contour_class::cont_both_junc;

            }
            else {
                // Otherwise the line has to be split. 
                for (l = 0; l <= counter; l++) {
                    begin = (l == 0) ? 0 : junc[i + l - 1].pos;
                    end = (l == counter) ? (tmp_cont.num - 1) : (l != counter && (i + l) > 0) ? junc[i + l].pos : 0;

                    if (end == begin && counter > 1)
                        continue;
#if 1
                    cont.emplace_back();
                    //till begin+num_pnt, since last index is exclusive in subList function
                    cont[num_cont].row = vector<float>(tmp_cont.row.begin(), tmp_cont.row.begin() + num_pnt);
                    cont[num_cont].col = vector<float>(tmp_cont.col.begin(), tmp_cont.col.begin() + num_pnt);
                    cont[num_cont].angle = vector<float>(tmp_cont.angle.begin(), tmp_cont.angle.begin() + num_pnt);
                    cont[num_cont].response = vector<float>(tmp_cont.response.begin(), tmp_cont.response.begin() + num_pnt);
                    cont[num_cont].num = num_pnt;
#endif
                    // Modify contour class. 
                    if (l == 0) {
                        if (tmp_cont.cont_class == contour_class::cont_start_junc ||
                            tmp_cont.cont_class == contour_class::cont_both_junc)
                            cont[(int)num_cont].cont_class = contour_class::cont_both_junc;
                        else
                            cont[(int)num_cont].cont_class = contour_class::cont_end_junc;
                    }
                    else if (l == counter) {
                        if (tmp_cont.cont_class == contour_class::cont_end_junc ||
                            tmp_cont.cont_class == contour_class::cont_both_junc)
                            cont[(int)num_cont].cont_class = contour_class::cont_both_junc;
                        else
                            cont[(int)num_cont].cont_class = contour_class::cont_start_junc;
                    }
                    /*      else {
                              cont[(int)num_cont].cont_class = contour_class::cont_both_junc;
                          }
                          num_cont++;
                      }
                      cont[j] = cont[--num_cont]; */
                }
            }

        }
    }



    // Finally, check whether all angles point to the right of the line. 
    for (i = 0; i < num_cont; i++)
    {
        tmp_cont = cont[i];
        if (tmp_cont.row.empty()) continue;

        num_pnt = tmp_cont.num;
        trow = tmp_cont.row;
        tcol = tmp_cont.col;
        tangle = tmp_cont.angle;

        // One point of the contour is enough to determine the orientation. 
        k = (int)((num_pnt - 1) / 2);
        if (k == 0 || trow.size() == 1 || k > (trow.size() - 2)) continue;

        // The next few lines are ok because lines have at least two points. 
        dx = trow[k + 1] - trow[k];
        dy = tcol[k + 1] - tcol[k];
        nx = std::cos(tangle[k]);
        ny = std::sin(tangle[k]);
        // If the angles point to the left of the line they have to be adapted.
        //  The orientation is determined by looking at the z-component of the
        //  cross-product of (dx,dy,0) and (nx,ny,0). 
        if (nx * dy - ny * dx < 0){
            for (j = 0; j < num_pnt; j++)
            {
                tangle[j] = (float)(tangle[j] + pi);
                if (tangle[j] >= 2 * pi)
					tangle[j] = (float)(tangle[j] - 2 * pi);
            }
        }
    }

    //Remove lines with number of points less than threshold
    m_results.resize(0);

    for (auto cc = 0; cc < cont.size(); cc++) {
        const contour& ct = cont[cc];
        auto len = ct.compute_length();
        if (nMinNumberOfNodes > 0 && len < nMinNumberOfNodes)
            continue;
        auto min_row = std::min_element(ct.row.begin(), ct.row.end());
        if (signbit(*min_row)) continue;
        auto min_col = std::min_element(ct.col.begin(), ct.col.end());
        if (signbit(*min_col)) continue;

        m_results.emplace_back(ct);
    }

    if (debug()) {
        auto actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        std::cout << "Computing Contours Time " << actual_wait_time << " microseconds" << std::endl;;
    }
}