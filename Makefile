SHELL = /bin/bash

#### Start of system configuration section. ####

srcdir = .
altsrcdir = ../samutil

prefix = $(HOME)/usr
bindir = $(prefix)/bin

C = gcc
CC = g++
INSTALL = /usr/bin/install -c
INSTALLDATA = /usr/bin/install -c -m 644
OBJDIR = obj
YEPLIBDIR = $(HOME)/cc/yeppp/library/binaries/x64-linux-sysv-default


GSLDEBUGLIB = /usr/lib
#GSLDEBUGLIB = /usr/lib/debug/usr/lib
#GSLDEBUGLIB= $(HOME)/usr/lib/

# This isn't set up to use 'yeppp/bla.h'.  Instead, the headers are
# prefixed with 'yep'
YEPDIR = $(HOME)/cc/yeppp/library/headers

# These are set so that we can say 'klib/bla.h' or 'htslib/bla.h'
SDGDIR = $(HOME)/cc/
HTSDIR = $(HOME)/cc/htslib
HTSLIBDIR = $(HTSDIR)

CPPFLAGS = -I. -I$(YEPDIR) -I$(SDGDIR) -I$(HTSDIR) $(EXTRA_CPPFLAGS)
OPT = -O0 -ggdb3
PROF = 
CXXFLAGS = $(OPT) $(PROF) -Wall -std=gnu++0x
CFLAGS = $(OPT) $(PROF) -Wall -std=gnu99
LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -lm -lgmp -lz -lpthread -lrt $(PROF)

DEPLIBS = -lgsl -lgslcblas -lm -lyeppp -lz -lpthread -lrt -lhts
#LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -llevmar -lm -lgmpxx -lgmp -llapack -lblas -lgfortran -lcblas -latlas

SOURCES = $(shell find $(srcdir) -name "*.cc")
CSOURCES = $(shell find $(srcdir) -name '*.c')

ALTSOURCES = $(shell find $(altsrcdir) -name "*.cc")

EXE = dep test_dirichlet
#EXE += quantile_test filter_sam_by_score

.PHONY : all
all : $(EXE)

# temporarily rename to dep_dev so as not to interfere with running binary.
dep : $(addprefix $(OBJDIR)/, dep.o common_tools.o dist.o				\
	locus_diff.o fasta.o locus_range.o wquantile.o sampling.o timer.o	\
	bam_sample_info.o binomial_est.o dir_cache.o						\
	dirichlet_points_gen.o dirichlet_diff_cache.o file_utils.o			\
	ordering.o locus.o bam_reader.o batch_pileup.o thread_queue.o		\
	virtual_bound.o gen_pair_comp.o geometry.o simplex.o				\
	chunk_strategy.o dep_pileup.o)
	$(CC) -L$(YEPLIBDIR) -L$(HTSLIBDIR)									\
	-Wl,-rpath,$(YEPLIBDIR),-rpath,$(GSLDEBUGLIB),-rpath,$(HTSLIBDIR)	\
	$(PROF) -o $@ $^ $(DEPLIBS)

test_distance : $(addprefix $(OBJDIR)/, test_distance.o spatial_search.o)
	$(C) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lrt -lm

test_likelihood : $(addprefix $(OBJDIR)/, test_likelihood.o likelihood.o)
	$(C) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

quantile_test : quantile_test.o sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

window_average : window_average.o histo.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lz -lrt

pileup_to_bindepth : obj/pileup_to_bindepth.o obj/bindepth.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lz -lrt

bindepth_to_pileup : obj/bindepth_to_pileup.o obj/bindepth.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lz -lrt

pileup_depth_stats : pileup_depth_stats.o bindepth.o histo.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrt

fasta_grep : fasta_grep.o bindepth.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_dirichlet : test_dirichlet.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

filter_sam_by_score : filter_sam_by_score.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

distance_wise : obj/distance_wise.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

-include $(patsubst ./%.cc,$(OBJDIR)/%.d,$(SOURCES))
-include $(patsubst ./%.c,$(OBJDIR)/%.d,$(CSOURCES))

define make-depend
$(CXX) -MM -MF $1 -MP -MT $2 $(CXXFLAGS) $(CPPFLAGS) $3
endef

define cmake-depend
$(C) -MM -MF $1 -MP -MT $2 $(CFLAGS) $(CPPFLAGS) $3
endef

$(OBJDIR)/%.o: %.cc
	$(call make-depend,$(subst .o,.d,$@),$@,$<)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.c
	$(call cmake-depend,$(subst .o,.d,$@),$@,$<)
	$(C) $(CFLAGS) $(CPPFLAGS) -c $< -o $@


.PHONY: install
install: $(EXE)
	$(INSTALL) $^ $(bindir)


.PHONY: clean
clean:
	rm -f *.o
