SMILEICXX     ?= mpic++
HDF5_ROOT_DIR ?=

PYTHONCONFIG ?= python-config
##################################################
DESCRIBE:=$(shell git describe 2>/dev/null || echo '??')
BRANCH:=$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '??')
COMMITDATE:=$(shell git show -s --pretty="%ci" 2>/dev/null || echo '??')

VERSION:="$(DESCRIBE)-$(BRANCH)"

#SuperLU_DIR = /home/wp/codes/science/SuperLU_4.3
LIB_Blas = $(SuperLU_DIR)/lib/blas_linux.a
LIB_SuperLU = $(SuperLU_DIR)/lib/libsuperlu_4.3.a
##LIBS = $(LIB_Blas) -L$(SuperLU_DIR)/lib -lsuperlu_4.3
LIBS = $(LIB_SuperLU) $(LIB_Blas)




CFLAGS +=  -D__VERSION=\"$(VERSION)\" -D__COMMITDATE=\""$(COMMITDATE)"\" -I${HDF5_ROOT_DIR}/include -std=c++0x
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5 -lz

CFLAGS += -I$(SuperLU_DIR)/SRC
LDFLAGS += $(LIBS) -lm

ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz -L/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/lib
endif

#add subdirs
DIRS = $(shell find * -type d)
#add include directives for subdirs
CFLAGS += $(DIRS:%=-I%)

EXEC = smilei

default: $(EXEC)

help:
	@echo 'Usage                 : make [clean] [openmp={gnu|intel}] [config=debug]'
	@echo 'Environment variables : SMILEICXX (c++ compiler) HDF5_ROOT_DIR (HDF5 dir with /include and /lib)'


# these are kept for backward compatibility and might be removed (see make help)
release:
	make

openmpgnu:
	make openmp=gnu

openmpintel:
	make openmp=intel

debug:
	make config=debug

scalasca:
	make config=scalasca


PYSCRIPTS := $(shell find Python -name \*.py)
PYHEADERS := $(PYSCRIPTS:.py=.pyh)

#collect all cpp files
SRCS := $(shell find . -name \*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)


CFLAGS+=$(shell $(PYTHONCONFIG) --includes)
LDFLAGS+=$(shell $(PYTHONCONFIG) --ldflags)


ifeq ($(config),debug)
	CFLAGS += -g -pg -Wall -D__DEBUG -O0# -shared-intel
else
	CFLAGS += -O3#  -xHost -ipo
endif

ifeq ($(openmp),gnu)
	CFLAGS += -fopenmp -D_OMP
else ifeq ($(openmp),intel)
	CFLAGS += -openmp -D_OMP
endif

ifeq ($(config),scalasca)
	CFLAGS += -openmp -D_OMP
    SMILEICXX = scalasca -instrument mpic++
endif

clean:
	rm -f $(OBJS) $(DEPS) $(PYHEADERS)

distclean: clean
	rm -f $(EXEC)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

# this generates a .h file containing a char[] with the python script in binary then
#you can just include this file to get the contents
%.pyh: %.py
	@ echo "Creating binary char for $< : $@"
	@ xxd -i $< > $@

%.d: %.cpp $(PYHEADERS)
	@ echo "Checking dependencies for $<"
# create and modify dependecy file .d to take into account the location subdir
	@ $(SMILEICXX) $(CFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(shell  dirname $<)/\1.d $(shell  dirname $<)/\1.o:@" > $@

%.o : %.cpp
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(PYHEADERS) $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

# DO NOT DELETE
