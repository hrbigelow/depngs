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
GSLDEBUGLIB = /usr/lib/debug/usr/lib
YEPHEADERS = $(HOME)/cc/yeppp/library/headers
CPPFLAGS = -I. -I.. -I$(YEPHEADERS)
OPT = -O0
PROF = 
CXXFLAGS = -ggdb3 $(OPT) $(PROF) -Wall -std=gnu++0x
CFLAGS = -ggdb3 $(OPT) $(PROF) -Wall -std=gnu99
LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -lm -lgmp -lz -lpthread -lrt

DEPLIBS = -lgsl -lgslcblas -lm -lyeppp -lz -lpthread
#LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -llevmar -lm -lgmpxx -lgmp -llapack -lblas -lgfortran -lcblas -latlas

SOURCES = $(shell find $(srcdir) -name "*.cc")
CSOURCES = $(shell find $(srcdir) -name '*.c')

ALTSOURCES = $(shell find $(altsrcdir) -name "*.cc")

EXE = dep test_dirichlet
#EXE += quantile_test filter_sam_by_score

.PHONY : all
all : $(EXE)

dep : $(addprefix $(OBJDIR)/, dep.o comp.o dict.o bqs.o bqs2jpd.o	\
	sampling.o tools.o nucleotide_stats.o pileup_tools.o			\
	metropolis_sampling.o usage_strings.o run_comp.o dist.o			\
	dist_worker.o comp_worker.o pug.o file_utils.o					\
	file_binary_search.o ordering.o locus.o range_line_reader.o		\
	thread_queue.o)
	$(CC) -L$(YEPLIBDIR) -L$(GSLDEBUGLIB) \
	-Wl,-rpath,$(YEPLIBDIR),-rpath,$(GSLDEBUGLIB) -o $@ $^ $(DEPLIBS)

test_distance : $(addprefix $(OBJDIR)/, test_distance.o spatial_search.o)
	$(C) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lrt -lm

test_likelihood : $(addprefix $(OBJDIR)/, test_likelihood.o likelihood.o)
	$(C) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

testopt : $(addprefix $(OBJDIR)/, pow_int.o testopt.o)
	$(C) $(CFLAGS) -L$(YEPLIBDIR) -o $@ $^ -lgsl -lgslcblas -lm -lyeppp



quantile_test : quantile_test.o sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

1dhist : 1dhist.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lz -lrt

window_average : window_average.o histo.o ../samutil/obj/file_utils.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lz -lrt

#strip_pileup : strip_pileup.o ../samutil/obj/file_utils.o
#	$(CXX) $(CXXFLAGS) -o $@ $^ -lz -lrt

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
