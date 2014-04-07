SHELL = /bin/bash

#### Start of system configuration section. ####

srcdir = .
altsrcdir = ../samutil

prefix = $(HOME)/usr
bindir = $(prefix)/bin

CC = g++
INSTALL = /usr/bin/install -c
INSTALLDATA = /usr/bin/install -c -m 644
OBJDIR = obj
CPPFLAGS = -I. -I..
OPT = -O0
PROF = 
CXXFLAGS = -ggdb3 $(OPT) $(PROF) -Wall -std=gnu++0x
LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -lm -lgmp -lz -lpthread -lrt

#LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -llevmar -lm -lgmpxx -lgmp -llapack -lblas -lgfortran -lcblas -latlas

SOURCES = $(shell find $(srcdir) -name "*.cc")
ALTSOURCES = $(shell find $(altsrcdir) -name "*.cc")

EXE = dep test_dirichlet
#EXE += quantile_test filter_sam_by_score

.PHONY : all
all : $(EXE)

dep : $(addprefix $(OBJDIR)/, dep.o comp.o mode.o comp_functor.o		\
	simp.o simc.o bqs.o bqslocus.o bqs2jpd.o metropolis.o sampling.o	\
	posterior.o tools.o transformation.o error_estimate.o				\
	pileup_tools.o stats_tools.o dirichlet.o slice_sampling.o			\
	hilbert.o simulation.o nucleotide_stats.o usage_strings.o			\
	run_comp_or_mode.o anomaly_tools.o dist.o dist_worker.o)			\
	../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

quantile_test : quantile_test.o sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)


test_dirichlet : test_dirichlet.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

filter_sam_by_score : filter_sam_by_score.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

distance_wise : obj/distance_wise.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

-include $(patsubst ./%.cc,$(OBJDIR)/%.d,$(SOURCES))

define make-depend
$(CXX) -MM -MF $3 -MP -MT $2 $(CXXFLAGS) $(CPPFLAGS) $1
endef


$(OBJDIR)/%.o: %.cc
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@


.PHONY: install
install: $(EXE)
	$(INSTALL) $^ $(bindir)


.PHONY: clean
clean:
	rm -f *.o
