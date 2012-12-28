SamplingTest : SamplingTest.o slice_sampling.o metropolis.o error_estimate.o hilbert.o tools.o sampling.o \
	integrands.o stats_tools.o simulation.o transformation.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)


HilbertTest : HilbertTest.o hilbert.o slice_sampling.o sampling.o integrands.o error_estimate.o tools.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

StepTest : StepTest.o hilbert.o slice_sampling.o tools.o sampling.o integrands.o error_estimate.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

MetropolisTest : MetropolisTest.o tools.o sampling.o integrands.o error_estimate.o metropolis.o stats_tools.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

