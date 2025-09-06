#pragma once
#include<cmath>

class BivariateQuadratic {
public:
    // Constructor to set grid parameters
    BivariateQuadratic(double dx_val = 1.0, double dy_val = 1.0, 
                      double xmin_val = 0.0, double xmax_val = 2.0,
                      double ymin_val = 0.0, double ymax_val = 2.0) 
        : dx(dx_val), dy(dy_val), xmin(xmin_val), xmax(xmax_val), 
          ymin(ymin_val), ymax(ymax_val) {}
    
    double interp_bilinear(double f3x3[3][3], double x, double y) {
        int i = (int)std::floor((x - xmin) / dx);
        int j = (int)std::floor((y - ymin) / dy);
        
        // Handle boundary case: when x = xmin + dx, we want i=0, s=1 (not i=1, s=0)
        if (x == xmin + dx) i = 0;
        if (y == ymin + dy) j = 0;
        
        double s = (x - xmin) / dx - i;
        double t = (y - ymin) / dy - j;

        return (1-t)*((1-s)*f3x3[i][j] + s*f3x3[i+1][j]) + t*((1-s)*f3x3[i][j+1] + s*f3x3[i+1][j+1]);
    }
    double interp_biquad(double f3x3[3][3], double x, double y) {
        double i = (int)std::floor((x-xmin)/dx);
        double j = (int)std::floor((y-ymin)/dy);
        
        // Handle boundary case: when x = xmin + dx, we want i=0, s=1 (not i=1, s=0)
        if (x == xmin + dx) i = 0;
        if (y == ymin + dy) j = 0;
        
        double s = (x - xmin) / dx - i;
        double t = (y - ymin) / dy - j;
        
        double Ls0 = (s - 1.0) * (s - 2.0) / 2.0;
        double Ls1 = s * (2.0 - s);
        double Ls2 = s * (s - 1.0) / 2.0;

        double intermediate0 = f3x3[0][0] * Ls0 + f3x3[1][0] * Ls1 + f3x3[2][0] * Ls2;
        double intermediate1 = f3x3[0][1] * Ls0 + f3x3[1][1] * Ls1 + f3x3[2][1] * Ls2;
        double intermediate2 = f3x3[0][2] * Ls0 + f3x3[1][2] * Ls1 + f3x3[2][2] * Ls2;

        double Lt0 = (t - 1.0) * (t - 2.0) / 2.0;
        double Lt1 = t * (2.0 - t);
        double Lt2 = t * (t - 1.0) / 2.0;
        return intermediate0 * Lt0 + intermediate1 * Lt1 + intermediate2 * Lt2;
    }
private:
    double dx, dy, xmin, xmax, ymin, ymax;
};
