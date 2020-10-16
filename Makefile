# Benchmark program for CSR SpMV
# Copyright (C) 2020 James D. Trotter
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Authors: James D. Trotter <james@simula.no>
#
# Benchmarking program for sparse matrix-vector multiplication (SpMV)
# with matrices in compressed sparse row (CSR) format.

csrspmv = csrspmv

all: $(csrspmv)
clean:
	rm -f $(csrspmv_c_objects) $(csrspmv)
.PHONY: all clean

CFLAGS += -g -Wall -iquote src

ifndef NO_OPENMP
CFLAGS += -fopenmp
endif

csrspmv_c_sources = \
	src/csr.c \
	src/main.c \
	src/matrix_market.c \
	src/parse.c \
	src/program_options.c \
	src/vector.c
csrspmv_c_headers = \
	src/csr.h \
	src/matrix_market.h \
	src/parse.h \
	src/program_options.h \
	src/vector.h
csrspmv_c_objects := $(foreach x,$(csrspmv_c_sources),$(x:.c=.o))
$(csrspmv_c_objects): %.o: %.c $(csrspmv_c_headers)
	$(CC) -c $(CFLAGS) $< -o $@
$(csrspmv): $(csrspmv_c_objects)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@
