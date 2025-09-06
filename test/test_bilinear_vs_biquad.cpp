#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include "../src/BivariateQuadratic.h"

// Test function that creates a grid from a known analytical function
void createGridFromFunction(double f3x3[3][3], std::function<double(double, double)> func, 
                          double xmin, double ymin, double dx, double dy) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double x = xmin + i * dx;
            double y = ymin + j * dy;
            f3x3[i][j] = func(x, y);
        }
    }
}

// Test that both methods exactly reproduce values at grid points within interpolation domain
bool testGridPointAccuracy(BivariateQuadratic& interp, double f3x3[3][3], 
                          std::function<double(double, double)> func,
                          double xmin, double ymin, double dx, double dy) {
    bool bilinear_passed = true;
    bool biquad_passed = true;
    
    // Test only grid points within the interpolation domain [xmin, xmin+dx) × [ymin, ymin+dy)
    for (int i = 0; i <= 1; i++) {  // Only i=0,1 (grid points at xmin, xmin+dx)
        for (int j = 0; j <= 1; j++) {  // Only j=0,1 (grid points at ymin, ymin+dy)
            double x = xmin + i * dx;
            double y = ymin + j * dy;
            double expected = func(x, y);
            
            double bilinear_result = interp.interp_bilinear(f3x3, x, y);
            double biquad_result = interp.interp_biquad(f3x3, x, y);
            
            double bilinear_error = std::abs(expected - bilinear_result);
            double biquad_error = std::abs(expected - biquad_result);
            
            if (bilinear_error > 1e-12) bilinear_passed = false;
            if (biquad_error > 1e-12) biquad_passed = false;
        }
    }
    
    return bilinear_passed && biquad_passed;
}

// Test polynomial reproduction properties
std::pair<double, double> testPolynomialReproduction(BivariateQuadratic& interp, 
                               std::function<double(double, double)> func,
                               double xmin, double ymin, double dx, double dy) {
    // Create grid from analytical function
    double f3x3[3][3];
    createGridFromFunction(f3x3, func, xmin, ymin, dx, dy);
    
    // Test at various points within the interpolation domain [xmin, xmin+dx) × [ymin, ymin+dy)
    // Use relative coordinates [0,1) and then scale to actual domain
    std::vector<std::pair<double, double>> relativeTestPoints = {
        {0.1, 0.1}, {0.3, 0.7}, {0.5, 0.5}, {0.8, 0.2}, {0.2, 0.8},
        {0.25, 0.25}, {0.75, 0.75}, {0.9, 0.1}, {0.1, 0.9}, {0.6, 0.4}
    };
    
    double max_bilinear_error = 0.0;
    double max_biquad_error = 0.0;
    
    for (const auto& relPoint : relativeTestPoints) {
        // Scale relative coordinates to actual domain
        double x = xmin + relPoint.first * dx;
        double y = ymin + relPoint.second * dy;
        // Test points should be within the interpolation domain [xmin, xmin+dx) × [ymin, ymin+dy)
        if (x >= xmin && x < xmin + dx && y >= ymin && y < ymin + dy) {
            double expected = func(x, y);
            double bilinear_result = interp.interp_bilinear(f3x3, x, y);
            double biquad_result = interp.interp_biquad(f3x3, x, y);
            
            double bilinear_error = std::abs(expected - bilinear_result);
            double biquad_error = std::abs(expected - biquad_result);
            
            max_bilinear_error = std::max(max_bilinear_error, bilinear_error);
            max_biquad_error = std::max(max_biquad_error, biquad_error);
        }
    }
    
    return {max_bilinear_error, max_biquad_error};
}

// Compare interpolation errors for non-polynomial functions
std::pair<double, double> testInterpolationAccuracy(BivariateQuadratic& interp,
                              std::function<double(double, double)> func,
                              double xmin, double ymin, double dx, double dy) {
    // Create grid from analytical function
    double f3x3[3][3];
    createGridFromFunction(f3x3, func, xmin, ymin, dx, dy);
    
    // Dense sampling for error analysis within the valid interpolation domain
    int n_samples = 21;  // 21x21 grid of test points
    double total_bilinear_error = 0.0;
    double total_biquad_error = 0.0;
    double max_bilinear_error = 0.0;
    double max_biquad_error = 0.0;
    int count = 0;
    
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_samples; j++) {
            // Sample within [xmin, xmin+dx) × [ymin, ymin+dy)
            double x = xmin + (double)i * dx / n_samples;  // Note: using n_samples instead of (n_samples-1)
            double y = ymin + (double)j * dy / n_samples;  // to stay within [xmin, xmin+dx)
            
            double expected = func(x, y);
            double bilinear_result = interp.interp_bilinear(f3x3, x, y);
            double biquad_result = interp.interp_biquad(f3x3, x, y);
            
            double bilinear_error = std::abs(expected - bilinear_result);
            double biquad_error = std::abs(expected - biquad_result);
            
            total_bilinear_error += bilinear_error;
            total_biquad_error += biquad_error;
            max_bilinear_error = std::max(max_bilinear_error, bilinear_error);
            max_biquad_error = std::max(max_biquad_error, biquad_error);
            count++;
        }
    }
    
    double avg_bilinear_error = total_bilinear_error / count;
    double avg_biquad_error = total_biquad_error / count;
    
    return {avg_bilinear_error, avg_biquad_error};
}

int main() {
    std::cout << "=================================================================\n";
    std::cout << "     BILINEAR vs BIQUADRATIC INTERPOLATION COMPARISON TESTS     \n";
    std::cout << "=================================================================\n";
    
    // Create interpolator with unit grid spacing
    double dx = 0.25, dy = 0.25;
    double xmin = 0.0, xmax = 0.5;
    double ymin = 0.0, ymax = 0.5;
    BivariateQuadratic interp(dx, dy, xmin, xmax, ymin, ymax);
    
    // Test functions of increasing complexity
    std::vector<std::pair<std::function<double(double, double)>, std::string>> testFunctions = {
        // Polynomials - should be exactly reproduced up to certain degrees
        {[](double, double) { return 1.0; }, "Constant (f=1)"},
        {[](double x, double) { return x; }, "Linear in x (f=x)"},
        {[](double, double y) { return y; }, "Linear in y (f=y)"},
        {[](double x, double y) { return x + y; }, "Bilinear (f=x+y)"},
        {[](double x, double y) { return x * y; }, "Product (f=xy)"},
        {[](double x, double) { return x*x; }, "Quadratic in x (f=x²)"},
        {[](double, double y) { return y*y; }, "Quadratic in y (f=y²)"},
        {[](double x, double y) { return x*x + y*y; }, "Quadratic sum (f=x²+y²)"},
        {[](double x, double y) { return x*x*y + x*y*y; }, "Mixed cubic (f=x²y+xy²)"},
        
        // Non-polynomial functions - show approximation differences
        {[](double x, double y) { return sin(M_PI*x/2) * sin(M_PI*y/2); }, "Trigonometric"},
        {[](double x, double y) { return exp(x*y); }, "Exponential"},
        {[](double x, double y) { return sqrt(x*x + y*y + 0.1); }, "Radial distance"},
        {[](double x, double y) { return 1.0 / (1.0 + x*x + y*y); }, "Rational"}
    };
    
    std::cout << "\nTest Results Summary:\n";
    std::cout << std::string(82, '-') << "\n";
    std::cout << std::left << std::setw(18) << "Function";
    std::cout << std::setw(10) << "Grid Pts";
    std::cout << std::setw(16) << "Poly Max Error";
    std::cout << std::setw(16) << "Interp Avg Err";
    std::cout << std::setw(10) << "Better";
    std::cout << "\n";
    std::cout << std::setw(18) << " ";
    std::cout << std::setw(10) << "(Pass?)";
    std::cout << std::setw(8) << "Bilin";
    std::cout << std::setw(8) << "Biquad";
    std::cout << std::setw(8) << "Bilin";
    std::cout << std::setw(8) << "Biquad";
    std::cout << std::setw(10) << "Method";
    std::cout << "\n";
    std::cout << std::string(82, '-') << "\n";
    
    for (size_t i = 0; i < testFunctions.size(); i++) {
        auto func = testFunctions[i].first;
        auto name = testFunctions[i].second;
        
        // Test grid point accuracy for all functions
        double f3x3[3][3];
        createGridFromFunction(f3x3, func, xmin, ymin, dx, dy);
        bool gridPass = testGridPointAccuracy(interp, f3x3, func, xmin, ymin, dx, dy);
        
        // Test polynomial reproduction
        std::pair<double, double> polyResults = testPolynomialReproduction(interp, func, xmin, ymin, dx, dy);
        double polyBilinErr = polyResults.first;
        double polyBiquadErr = polyResults.second;
        
        // Test interpolation accuracy
        std::pair<double, double> interpResults = testInterpolationAccuracy(interp, func, xmin, ymin, dx, dy);
        double interpBilinErr = interpResults.first;
        double interpBiquadErr = interpResults.second;
        
        // Print results in table format
        std::cout << std::left << std::setw(18) << name.substr(0, 17);
        std::cout << std::setw(10) << (gridPass ? "PASS" : "FAIL");
        std::cout << std::scientific << std::setprecision(1);
        std::cout << std::setw(8) << polyBilinErr;
        std::cout << std::setw(8) << polyBiquadErr;
        std::cout << std::setw(8) << interpBilinErr;
        std::cout << std::setw(8) << interpBiquadErr;
        std::cout << std::setw(10) << (interpBiquadErr < interpBilinErr ? "Biquad" : "Bilinear");
        std::cout << "\n";
    }
    
    return 0;
}
