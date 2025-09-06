CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
SRCDIR = src
TESTDIR = test
BUILDDIR = build

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile the tests
$(BUILDDIR)/test_version2_correctness: $(TESTDIR)/test_version2_correctness.cpp $(SRCDIR)/BivariateQuadratic.h | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I. $(TESTDIR)/test_version2_correctness.cpp -o $(BUILDDIR)/test_version2_correctness

$(BUILDDIR)/test_bilinear_correctness: $(TESTDIR)/test_bilinear_correctness.cpp $(SRCDIR)/BivariateQuadratic.h | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I. $(TESTDIR)/test_bilinear_correctness.cpp -o $(BUILDDIR)/test_bilinear_correctness

$(BUILDDIR)/test_bilinear_vs_biquad: $(TESTDIR)/test_bilinear_vs_biquad.cpp $(SRCDIR)/BivariateQuadratic.h | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I. $(TESTDIR)/test_bilinear_vs_biquad.cpp -o $(BUILDDIR)/test_bilinear_vs_biquad

# Run the version2 correctness test  
test_biquad: $(BUILDDIR)/test_version2_correctness
	./$(BUILDDIR)/test_version2_correctness

# Run the bilinear correctness test
test_bilinear: $(BUILDDIR)/test_bilinear_correctness
	./$(BUILDDIR)/test_bilinear_correctness

# Run the bilinear vs biquadratic comparison test
test_comparison: $(BUILDDIR)/test_bilinear_vs_biquad
	./$(BUILDDIR)/test_bilinear_vs_biquad

# Run all tests
test: test_biquad test_bilinear test_comparison

# Clean build files
clean:
	rm -rf $(BUILDDIR)
	rm -rf *.dat

.PHONY: test test_biquad test_bilinear test_comparison clean
