#! /usr/bin/python3

import libstreams
import mpi4py

print(dir(libstreams))

libstreams.wrap_startmpi()
libstreams.wrap_setup()

# now the main solver loop
libstreams.wrap_init_solver()

