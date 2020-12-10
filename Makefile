CXXFLAGS     = -std=c++14 -g -fopenmp -Ofast -mavx2 -mfma
OBJDIR       = obj
DEPDIR       = $(OBJDIR)/.deps
# Flags which, when added to gcc/g++, will auto-generate dependency files
DEPFLAGS     = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

# Function which takes a list of words and returns a list of unique words in that list
# https://stackoverflow.com/questions/16144115/makefile-remove-duplicate-words-without-sorting
uniq         = $(if $1,$(firstword $1) $(call uniq,$(filter-out $(firstword $1),$1)))

# Source files - add more to auto-compile into .o files
SOURCES      = Common/image.cpp Common/fft.cpp Common/mask.cpp Experiment1/main.cpp Experiment2/main.cpp Experiment3/main.cpp
# Executable targets - add more to auto-make in default 'all' target
EXEC         = Experiment1/remove-noise Experiment2/frequency-filter Experiment3/homomorphic
# Targets required for the homework, spearated by experiment
REQUIRED_1   = out/boy.pgm
REQUIRED_2   = 
REQUIRED_3   = out/girl-0.5-1.5-homomorphic.pgm
REQUIRED_OUT = $(REQUIRED_1) $(REQUIRED_2) $(REQUIRED_3)

SOURCEDIRS   = $(call uniq, $(dir $(SOURCES)))
OBJDIRS      = $(addprefix $(OBJDIR)/, $(SOURCEDIRS))
DEPDIRS      = $(addprefix $(DEPDIR)/, $(SOURCEDIRS))
DEPFILES     = $(SOURCES:%.cpp=$(DEPDIR)/%.d)

.PHONY: all clean report

# By default, make all executable targets and the images required for the homework
all: $(EXEC) $(REQUIRED_OUT)

# Executable Targets
Experiment1/remove-noise: $(OBJDIR)/Experiment1/main.o $(OBJDIR)/Common/image.o $(OBJDIR)/Common/fft.o
	$(CXX) $(CXXFLAGS) $^ -o $@

Experiment2/frequency-filter: $(OBJDIR)/Experiment2/main.o $(OBJDIR)/Common/image.o $(OBJDIR)/Common/fft.o $(OBJDIR)/Common/mask.o
	$(CXX) $(CXXFLAGS) $^ -o $@

Experiment3/homomorphic: $(OBJDIR)/Experiment3/main.o $(OBJDIR)/Common/image.o $(OBJDIR)/Common/fft.o $(OBJDIR)/Common/mask.o
	$(CXX) $(CXXFLAGS) $^ -o $@

### Experiment 1 Outputs ###
out/boy.pgm out/boy_smoothed.pgm out/boy_noise.pgm out/boy_spectrum.png: Experiment1/remove-noise | out
	Experiment1/remove-noise Images/boy_noisy.pgm out/boy.pgm

### Experiment 2 Outputs ###
out/lenna_filtered_frequency.pgm out/lenna_filtered_spatial.pgm out/sobel_imag.dat out/sobel_spectrum.pgm: Experiment2/frequency-filter | out
	Experiment2/frequency-filter Images/lenna.pgm out/lenna_filtered.pgm

### Experiment 3 Outputs ###
.SECONDEXPANSION:
out/%-homomorphic.pgm: Experiment3/homomorphic Images/$$(word 1,$$(subst -, ,$$*)).pgm | out
	Experiment3/homomorphic Images/$(word 1,$(subst -, ,$*)).pgm -gl $(word 2,$(subst -, ,$*)) -gh $(word 3,$(subst -, ,$*)) $@

out/filter-%.pdf: Experiment3/filter.plt | out
	gnuplot -e "c=$*" -e "outfile='$@'" Experiment3/filter.plt

# Figures needed for the report
report: Images/boy_noisy.png out/boy.png out/boy_smoothed_7.png out/boy_smoothed_15.png out/boy_noise.png out/boy_spectrum.png
report: out/lenna_filtered_spatial.png out/lenna_filtered_frequency.png
report: out/girl-0.5-1.5-homomorphic.png out/girl-0-2-homomorphic.png out/girl-1-1-homomorphic.png Images/girl.png
report: out/filter-1.0.pdf out/filter-0.5.pdf out/filter-3.0.pdf

clean:
	rm -rf $(OBJDIR)
	rm -f $(EXEC)
	rm -rf out
	rm -f Images/*.png

# Generate .png images from .pgm images. Needed for report, since pdfLaTeX doesn't support .pgm images
%.png: %.pgm
	pnmtopng $< > $@

# Auto-Build .cpp files into .o
$(OBJDIR)/%.o: %.cpp
$(OBJDIR)/%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIRS) $(OBJDIRS)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) -c $< -o $@

# Make generated directories
$(DEPDIRS) $(OBJDIRS) out: ; @mkdir -p $@
$(DEPFILES):
include $(wildcard $(DEPFILES))
