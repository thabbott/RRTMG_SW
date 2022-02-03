# Set search path and source file list
DIRS := $(shell cat Filepath)
VPATH := $(foreach dir,$(DIRS),$(wildcard $(dir)))
SOURCES := $(shell cat Srcfiles)
OBJS := $(addsuffix .o, $(basename $(SOURCES)))
.SUFFIXES: .f90 .o

all: $(BASE)/rrtmg_sw_wrapper.f90 $(OBJS)
	$(F2PY) -I. $(OBJS) -m rrtmg_sw_wrapper $(BASE)/rrtmg_sw_wrapper.f90

.f90.o:
	$(FC) $(FFLAGS) $<

include Depends
