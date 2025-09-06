#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
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

// Test that bilinear exactly reproduces the analytical function at grid points within interpolation domain
bool testGridPointReproduction(BivariateQuadratic& interp, double f3x3[3][3], 
                               std::function<double(double, double)> func,
                               double xmin, double ymin, double dx, double dy,
                               const std::string& funcName) {
    bool allPassed = true;
    int failCount = 0;
    
    // Only test grid points within the interpolation domain
    for (int i = 0; i <= 1; i++) {  // i = 0, 1 (points at xmin, xmin+dx)
        for (int j = 0; j <= 1; j++) {  // j = 0, 1 (points at ymin, ymin+dy)
            double x = xmin + i * dx;
            double y = ymin + j * dy;
            double expected = func(x, y);
            double interpolated = interp.interp_bilinear(f3x3, x, y);
            double diff = std::abs(expected - interpolated);
            
            if (diff >= 1e-12) {
                if (failCount == 0) {
                    std::cout << "Grid point reproduction FAILED for " << funcName << ":\n";
                }
                std::cout << "  Point (" << x << "," << y << "): expected=" << expected 
                         << ", got=" << interpolated << ", error=" << diff << "\n";
                allPassed = false;
                failCount++;
            }
        }
    }
    return allPassed;
}

// Test polynomial reproduction property (bilinear should exactly reproduce polynomials up to degree 1)
bool testPolynomialReproduction(BivariateQuadratic& interp, 
                               std::function<double(double, double)> func,
                               const std::string& funcName,
                               double xmin, double ymin, double dx, double dy) {
    // Create grid from analytical function
    double f3x3[3][3];
    createGridFromFunction(f3x3, func, xmin, ymin, dx, dy);
    
    // Test at various points within the interpolation domain [xmin, xmin+dx] × [ymin, ymin+dy]
    std::vector<std::pair<double, double>> testPoints = {
        {0.1, 0.1}, {0.3, 0.7}, {0.5, 0.5}, {0.8, 0.2}, {0.2, 0.8},
        {0.15, 0.15}, {0.7, 0.9}, {0.25, 0.75}, {0.9, 0.1}, {0.6, 0.4}
    };
    
    bool allPassed = true;
    double maxError = 0.0;
    int failCount = 0;
    
    for (const auto& point : testPoints) {
        double x = point.first, y = point.second;
        // Test points should be within the interpolation domain [xmin, xmin+dx] × [ymin, ymin+dy]
        if (x >= xmin && x <= xmin + dx && y >= ymin && y <= ymin + dy) {
            double expected = func(x, y);
            double interpolated = interp.interp_bilinear(f3x3, x, y);
            double diff = std::abs(expected - interpolated);
            maxError = std::max(maxError, diff);
            
            if (diff >= 1e-10) {
                if (failCount == 0) {
                    std::cout << "Polynomial reproduction FAILED for " << funcName << ":\n";
                }
                std::cout << "  Point (" << std::setw(4) << x << "," << std::setw(4) << y << "): "
                         << "expected=" << std::setw(10) << std::setprecision(6) << expected 
                         << ", got=" << std::setw(10) << interpolated 
                         << ", error=" << std::scientific << std::setprecision(2) << diff << "\n";
                allPassed = false;
                failCount++;
            }
        }
    }
    
    if (!allPassed) {
        std::cout << std::fixed << "  Maximum error: " << std::scientific << maxError << "\n";
    }
    return allPassed;
}

// Test bilinear basis function properties
bool testBilinearBasisProperties() {
    // Bilinear basis functions: (1-s)(1-t), s(1-t), (1-s)t, st
    // where s,t ∈ [0,1] are normalized coordinates
    
    // Test that basis functions sum to 1 (partition of unity)
    std::vector<std::pair<double, double>> testPoints = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0},  // corners
        {0.5, 0.5}, {0.25, 0.75}, {0.8, 0.2}, {0.1, 0.9}  // interior points
    };
    
    bool passed = true;
    int failCount = 0;
    
    for (const auto& point : testPoints) {
        double s = point.first, t = point.second;
        
        // Bilinear basis functions
        double phi00 = (1.0 - s) * (1.0 - t);
        double phi10 = s * (1.0 - t);
        double phi01 = (1.0 - s) * t;
        double phi11 = s * t;
        
        double sum = phi00 + phi10 + phi01 + phi11;
        double diff = std::abs(sum - 1.0);
        
        if (diff >= 1e-15) {
            if (failCount == 0) {
                std::cout << "Bilinear basis partition of unity FAILED:\n";
            }
            std::cout << "  s=" << s << ", t=" << t << ": sum = " << sum 
                     << ", error from 1 = " << std::scientific << diff << "\n";
            passed = false;
            failCount++;
        }
    }
    
    // Test that basis functions have correct values at corners
    std::vector<std::vector<double>> corners = {{0,0}, {1,0}, {0,1}, {1,1}};
    std::vector<std::vector<double>> expected_values = {
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}
    };
    
    for (int corner = 0; corner < 4; corner++) {
        double s = corners[corner][0];
        double t = corners[corner][1];
        
        std::vector<double> phi = {
            (1.0 - s) * (1.0 - t),  // phi00
            s * (1.0 - t),          // phi10
            (1.0 - s) * t,          // phi01
            s * t                   // phi11
        };
        
        for (int i = 0; i < 4; i++) {
            double diff = std::abs(phi[i] - expected_values[corner][i]);
            if (diff >= 1e-15) {
                if (failCount == 0) {
                    std::cout << "Bilinear basis corner values FAILED:\n";
                }
                std::cout << "  Corner (" << s << "," << t << "), basis " << i 
                         << ": got=" << phi[i] << ", expected=" << expected_values[corner][i] 
                         << ", error=" << diff << "\n";
                passed = false;
                failCount++;
            }
        }
    }
    
    return passed;
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "=== Bilinear Interpolation Correctness Tests ===\n\n";
    
    // Test parameters
    double dx = 1.0, dy = 1.0;
    double xmin = 0.0, xmax = 2.0;
    double ymin = 0.0, ymax = 2.0;
    
    std::cout << "Interpolation domain: [" << xmin << ", " << xmin+dx << "] × [" << ymin << ", " << ymin+dy << "]\n\n";
    
    BivariateQuadratic interpolator(dx, dy, xmin, xmax, ymin, ymax);
    
    int testsPassed = 0;
    int totalTests = 0;
    
    // Test 1: Bilinear basis function properties
    bool basisTest = testBilinearBasisProperties();
    totalTests++;
    if (basisTest) testsPassed++;
    else std::cout << "Test 1: Bilinear Basis Properties - FAILED\n\n";
    
    // Test 2: Constant function (degree 0)
    auto constant = [](double, double) { return 42.0; };
    bool constTest = testPolynomialReproduction(interpolator, constant, "f(x,y) = 42", xmin, ymin, dx, dy);
    totalTests++;
    if (constTest) testsPassed++;
    else std::cout << "Test 2: Constant Function - FAILED\n\n";
    
    // Test 3: Linear functions (degree 1) - these should be exactly reproduced
    auto linear1 = [](double x, double y) { return 3.0 * x + 2.0 * y + 1.0; };
    bool linear1Test = testPolynomialReproduction(interpolator, linear1, "f(x,y) = 3x + 2y + 1", xmin, ymin, dx, dy);
    totalTests++;
    if (linear1Test) testsPassed++;
    else std::cout << "Test 3a: Linear Function 1 - FAILED\n\n";
    
    auto linear2 = [](double x, double y) { return -x + 4.0 * y - 2.0; };
    bool linear2Test = testPolynomialReproduction(interpolator, linear2, "f(x,y) = -x + 4y - 2", xmin, ymin, dx, dy);
    totalTests++;
    if (linear2Test) testsPassed++;
    else std::cout << "Test 3b: Linear Function 2 - FAILED\n\n";
    
    auto linearX = [](double x, double) { return 5.0 * x; };
    bool linearXTest = testPolynomialReproduction(interpolator, linearX, "f(x,y) = 5x", xmin, ymin, dx, dy);
    totalTests++;
    if (linearXTest) testsPassed++;
    else std::cout << "Test 3c: Linear X only - FAILED\n\n";
    
    auto linearY = [](double, double y) { return 3.0 * y; };
    bool linearYTest = testPolynomialReproduction(interpolator, linearY, "f(x,y) = 3y", xmin, ymin, dx, dy);
    totalTests++;
    if (linearYTest) testsPassed++;
    else std::cout << "Test 3d: Linear Y only - FAILED\n\n";
    
    // Test 4: Grid point reproduction for various functions
    double f3x3[3][3];
    
    createGridFromFunction(f3x3, constant, xmin, ymin, dx, dy);
    bool gridConst = testGridPointReproduction(interpolator, f3x3, constant, xmin, ymin, dx, dy, "constant");
    totalTests++;
    if (gridConst) testsPassed++;
    
    createGridFromFunction(f3x3, linear1, xmin, ymin, dx, dy);
    bool gridLinear = testGridPointReproduction(interpolator, f3x3, linear1, xmin, ymin, dx, dy, "linear");
    totalTests++;
    if (gridLinear) testsPassed++;
    
    // Test 5: Quadratic functions (should NOT be exactly reproduced - this tests approximation quality)
    auto quadX = [](double x, double) { return x * x; };
    
    double f3x3_quad[3][3];
    createGridFromFunction(f3x3_quad, quadX, xmin, ymin, dx, dy);
    
    std::vector<std::pair<double, double>> quadTestPoints = {
        {0.5, 0.5}, {0.25, 0.75}, {0.75, 0.25}, {0.1, 0.9}
    };
    
    double maxQuadError = 0.0;
    bool hasQuadErrors = false;
    std::cout << "Test 5: Quadratic Approximation Quality (x²):\n";
    for (const auto& point : quadTestPoints) {
        double x = point.first, y = point.second;
        if (x >= xmin && x <= xmin + dx && y >= ymin && y <= ymin + dy) {
            double expected = quadX(x, y);
            double interpolated = interpolator.interp_bilinear(f3x3_quad, x, y);
            double diff = std::abs(expected - interpolated);
            maxQuadError = std::max(maxQuadError, diff);
            std::cout << "  Point (" << x << "," << y << "): expected=" << expected 
                     << ", got=" << interpolated << ", error=" << diff << "\n";
            if (diff > 0.01) hasQuadErrors = true;  // Expect significant errors for quadratic
        }
    }
    std::cout << "  Maximum quadratic error: " << maxQuadError << "\n";
    
    if (!hasQuadErrors) {
        std::cout << "WARNING: Quadratic approximation errors unexpectedly small\n";
    }
    std::cout << "\n";
    
    // Test 6: Bilinear term xy (should NOT be exactly reproduced)
    auto bilinearTerm = [](double x, double y) { return x * y; };
    
    double f3x3_xy[3][3];
    createGridFromFunction(f3x3_xy, bilinearTerm, xmin, ymin, dx, dy);
    
    double maxXYError = 0.0;
    std::cout << "Test 6: Bilinear Term Approximation Quality (xy):\n";
    for (const auto& point : quadTestPoints) {
        double x = point.first, y = point.second;
        if (x >= xmin && x <= xmin + dx && y >= ymin && y <= ymin + dy) {
            double expected = bilinearTerm(x, y);
            double interpolated = interpolator.interp_bilinear(f3x3_xy, x, y);
            double diff = std::abs(expected - interpolated);
            maxXYError = std::max(maxXYError, diff);
            std::cout << "  Point (" << x << "," << y << "): expected=" << expected 
                     << ", got=" << interpolated << ", error=" << diff << "\n";
        }
    }
    std::cout << "  Maximum xy error: " << maxXYError << "\n\n";
    
    // Final summary
    std::cout << "=== SUMMARY ===\n";
    std::cout << "Tests passed: " << testsPassed << "/" << totalTests;
    
    if (testsPassed == totalTests) {
        std::cout << " - ALL TESTS PASSED! ✓\n";
        std::cout << "Bilinear interpolation correctly reproduces polynomials up to degree 1.\n";
        std::cout << "Quadratic approximation quality: max error = " << std::max(maxQuadError, maxXYError) << "\n";
        return 0;
    } else {
        std::cout << " - " << (totalTests - testsPassed) << " FAILED ✗\n";
        return 1;
    }
}
