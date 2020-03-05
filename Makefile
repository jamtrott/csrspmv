##
## Example program for sparse matrix-vector multiplication with the
## compressed sparse row (CSR) format.
##
programs = csr_spmv

clean-programs = $(programs:%=%-clean)

all: $(libs) $(programs)

clean: $(clean-programs)

.PHONY: all clean

##
## Configuration
##

INCLUDES = -iquote src
CFLAGS += -g -Wall
LDFLAGS += -lm

ifndef NO_OPENMP
CFLAGS += -fopenmp
endif

ifndef NO_LIBPFM
CFLAGS += -DHAVE_LIBPFM
LDFLAGS += -lpfm
endif


# An example program for computing a sparse matrix-vector
# multiplication with a matrix in CSR format.
csr_spmv = csr_spmv
csr_spmv_c_sources = \
	csr_spmv.c \
	mmio.c \
	perf_events.c \
	perf_session.c \
	sample_statistics.c
csr_spmv_c_headers = \
	mmio.h \
	perf_events.h \
	perf_session.h \
	sample_statistics.h
csr_spmv_c_objects := $(foreach x,$(csr_spmv_c_sources),$(x:.c=.o))
$(csr_spmv_c_objects): %.o: %.c $(csr_spmv_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(csr_spmv): $(csr_spmv_c_objects)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) -o $@
csr_spmv-clean:
	rm -f $(csr_spmv_c_objects) $(csr_spmv)
