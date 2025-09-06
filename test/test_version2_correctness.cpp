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

// Test that version2 exactly reproduces the analytical function at grid points within interpolation domain
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
            double interpolated = interp.interp_biquad(f3x3, x, y);
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

// Test polynomial reproduction property
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
            double interpolated = interp.interp_biquad(f3x3, x, y);
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

// Test Lagrange polynomial properties
bool testLagrangeProperties() {
    // Test that Lagrange polynomials satisfy L_i(j) = δ_ij
    auto L0 = [](double s) { return (s - 1.0) * (s - 2.0) / 2.0; };
    auto L1 = [](double s) { return s * (2.0 - s); };
    auto L2 = [](double s) { return s * (s - 1.0) / 2.0; };
    
    std::vector<double> testPoints = {0.0, 1.0, 2.0};
    bool passed = true;
    int failCount = 0;
    
    // Test Kronecker delta property L_i(j) = δ_ij
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double s = testPoints[j];
            double L_val;
            if (i == 0) L_val = L0(s);
            else if (i == 1) L_val = L1(s);
            else L_val = L2(s);
            
            double expected = (i == j) ? 1.0 : 0.0;
            double diff = std::abs(L_val - expected);
            
            if (diff >= 1e-15) {
                if (failCount == 0) {
                    std::cout << "Lagrange polynomial Kronecker delta property FAILED:\n";
                }
                std::cout << "  L" << i << "(" << s << ") = " << L_val 
                         << ", expected = " << expected << ", error = " << diff << "\n";
                passed = false;
                failCount++;
            }
        }
    }
    
    // Test partition of unity: L0 + L1 + L2 = 1
    std::vector<double> arbitrary_points = {-0.5, 0.3, 0.7, 1.2, 1.8, 2.5};
    
    for (double s : arbitrary_points) {
        double sum = L0(s) + L1(s) + L2(s);
        double diff = std::abs(sum - 1.0);
        
        if (diff >= 1e-15) {
            if (failCount == 0) {
                std::cout << "Lagrange polynomial partition of unity FAILED:\n";
            }
            std::cout << "  s=" << s << ": L0+L1+L2 = " << sum 
                     << ", error from 1 = " << std::scientific << diff << "\n";
            passed = false;
            failCount++;
        }
    }
    
    return passed;
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "=== Biquadratic Interpolation Correctness Tests ===\n\n";
    
    // Test parameters
    double dx = 1.0, dy = 1.0;
    double xmin = 0.0, xmax = 2.0;
    double ymin = 0.0, ymax = 2.0;
    
    std::cout << "Interpolation domain: [" << xmin << ", " << xmin+dx << "] × [" << ymin << ", " << ymin+dy << "]\n\n";
    
    BivariateQuadratic interpolator(dx, dy, xmin, xmax, ymin, ymax);
    
    int testsPassed = 0;
    int totalTests = 0;
    
    // Test 1: Lagrange polynomial properties
    bool lagrangeTest = testLagrangeProperties();
    totalTests++;
    if (lagrangeTest) testsPassed++;
    else std::cout << "Test 1: Lagrange Properties - FAILED\n\n";
    
    // Test 2: Constant function (degree 0)
    auto constant = [](double, double) { return 42.0; };
    bool constTest = testPolynomialReproduction(interpolator, constant, "f(x,y) = 42", xmin, ymin, dx, dy);
    totalTests++;
    if (constTest) testsPassed++;
    else std::cout << "Test 2: Constant Function - FAILED\n\n";
    
    // Test 3: Linear functions (degree 1)
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
    
    // Test 4: Pure quadratic terms (degree 2)
    auto quadX = [](double x, double) { return x * x; };
    bool quadXTest = testPolynomialReproduction(interpolator, quadX, "f(x,y) = x²", xmin, ymin, dx, dy);
    totalTests++;
    if (quadXTest) testsPassed++;
    else std::cout << "Test 4a: Quadratic x² - FAILED\n\n";
    
    auto quadY = [](double, double y) { return y * y; };
    bool quadYTest = testPolynomialReproduction(interpolator, quadY, "f(x,y) = y²", xmin, ymin, dx, dy);
    totalTests++;
    if (quadYTest) testsPassed++;
    else std::cout << "Test 4b: Quadratic y² - FAILED\n\n";
    
    auto quadXY = [](double x, double y) { return x * y; };
    bool quadXYTest = testPolynomialReproduction(interpolator, quadXY, "f(x,y) = xy", xmin, ymin, dx, dy);
    totalTests++;
    if (quadXYTest) testsPassed++;
    else std::cout << "Test 4c: Quadratic xy - FAILED\n\n";
    
    // Test 5: Full quadratic polynomial (degree 2)
    auto fullQuad = [](double x, double y) { 
        return 2.0*x*x + 3.0*y*y + 4.0*x*y + 5.0*x + 6.0*y + 7.0; 
    };
    bool fullQuadTest = testPolynomialReproduction(interpolator, fullQuad, 
        "f(x,y) = 2x² + 3y² + 4xy + 5x + 6y + 7", xmin, ymin, dx, dy);
    totalTests++;
    if (fullQuadTest) testsPassed++;
    else std::cout << "Test 5: Full Quadratic - FAILED\n\n";
    
    // Test 6: Grid point reproduction for various functions
    double f3x3[3][3];
    
    createGridFromFunction(f3x3, constant, xmin, ymin, dx, dy);
    bool gridConst = testGridPointReproduction(interpolator, f3x3, constant, xmin, ymin, dx, dy, "constant");
    totalTests++;
    if (gridConst) testsPassed++;
    
    createGridFromFunction(f3x3, linear1, xmin, ymin, dx, dy);
    bool gridLinear = testGridPointReproduction(interpolator, f3x3, linear1, xmin, ymin, dx, dy, "linear");
    totalTests++;
    if (gridLinear) testsPassed++;
    
    createGridFromFunction(f3x3, fullQuad, xmin, ymin, dx, dy);
    bool gridQuad = testGridPointReproduction(interpolator, f3x3, fullQuad, xmin, ymin, dx, dy, "quadratic");
    totalTests++;
    if (gridQuad) testsPassed++;
    
    // Test 7: Higher degree polynomial (should not be exact)
    auto cubic = [](double x, double y) { return x*x*x + y*y*y + x*x*y + x*y*y; };
    
    double f3x3_cubic[3][3];
    createGridFromFunction(f3x3_cubic, cubic, xmin, ymin, dx, dy);
    
    std::vector<std::pair<double, double>> cubicTestPoints = {
        {0.5, 0.5}, {0.15, 0.15}, {0.25, 0.75}, {0.75, 0.25}
    };
    
    double maxCubicError = 0.0;
    bool hasCubicErrors = false;
    for (const auto& point : cubicTestPoints) {
        double x = point.first, y = point.second;
        double expected = cubic(x, y);
        double interpolated = interpolator.interp_biquad(f3x3_cubic, x, y);
        double diff = std::abs(expected - interpolated);
        maxCubicError = std::max(maxCubicError, diff);
        if (diff > 0.1) hasCubicErrors = true;
    }
    
    if (!hasCubicErrors) {
        std::cout << "WARNING: Cubic approximation errors unexpectedly small (max: " << maxCubicError << ")\n\n";
    }
    
    // Final summary
    std::cout << "=== SUMMARY ===\n";
    std::cout << "Tests passed: " << testsPassed << "/" << totalTests;
    
    if (testsPassed == totalTests) {
        std::cout << " - ALL TESTS PASSED! ✓\n";
        std::cout << "Biquadratic Interpolation correctly implements tensor-product quadratic Lagrange interpolation.\n";
        return 0;
    } else {
        std::cout << " - " << (totalTests - testsPassed) << " FAILED ✗\n";
        return 1;
    }
}
