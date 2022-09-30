#! /usr/bin/python3

import libstreams
import mpi4py

print(dir(libstreams))

libstreams.wrap_startmpi()
libstreams.wrap_setup()

# now the main solver loop
libstreams.wrap_init_solver()


steps = 1

for i in range(steps):
    print(f"step {i}")
    libstreams.wrap_step_solver()

libstreams.wrap_finalize_solver()

print("finalizing solver")
libstreams.wrap_finalize()
