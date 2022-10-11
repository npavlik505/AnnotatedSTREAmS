#! /usr/bin/python3
print("RUNNING PYTHON VERSION")

import libstreams
libstreams.wrap_startmpi()

import mpi4py
import json


import io_utils
from config import Config

with open("/input/input.json", "r") as f:
    json_data = json.load(f)
    config = Config.from_json(json_data)

print(dir(libstreams))

libstreams.wrap_setup()

# initialize files
flowfields = io_utils.IoFile("/distribute_save/flowfields.h5")

# initialize datasets from files
grid_shape = [config.grid.nx, config.grid.ny, config.grid.nz]
velocity = flowfields.create_dataset("velocity", 1, [5, *grid_shape], 1)

# now the main solver loop
libstreams.wrap_init_solver()

libstreams.wrap_copy_gpu_to_cpu()

for i in range(config.temporal.num_iter):
    #print(f"step {i}")
    libstreams.wrap_step_solver()

# copy to CPU
libstreams.wrap_copy_gpu_to_cpu()
velocity.write_array(config.slice_flowfield_array(libstreams.mod_streams.w))


libstreams.wrap_finalize_solver()

print("finalizing solver")
libstreams.wrap_finalize()
