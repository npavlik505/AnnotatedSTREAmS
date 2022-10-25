#! /usr/bin/python3
print("RUNNING PYTHON VERSION")

import libstreams as streams 
streams.wrap_startmpi()

from mpi4py import MPI
import json

import io_utils
import numpy as np
from config import Config

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

is_master = rank == 0;

with open("/input/input.json", "r") as f:
    json_data = json.load(f)
    config = Config.from_json(json_data)


def hprint(*args):
    if is_master:
        print(*args)

span_average = np.zeros([5, config.nx_mpi(), config.ny_mpi()], dtype=np.float64)
temp_field = np.zeros((config.nx_mpi(), config.ny_mpi(), config.nz_mpi()), dtype=np.float64)

def calculate_span_averages(span_average, temp_field, streams_data):
    # for rhou, rhov, rhow, and rhoE
    for data_idx in [1,2,3,4]:
        temp_field[:] = np.divide(streams_data[data_idx, :, :, :], streams_data[0, :, :, :], out = temp_field[:])

        span_average[data_idx, :, :] = np.sum(temp_field, axis=2, out=span_average[data_idx, :, :])
        span_average[data_idx, :, :] /= config.grid.nz

    span_average[0, :, :] = np.sum(streams_data[0, :, :, :], axis=2, out=span_average[0, :, :])
    span_average[0, :, :] /= config.grid.nz


streams.wrap_setup()

# initialize files
flowfields = io_utils.IoFile("/distribute_save/flowfields.h5")

# initialize datasets from files
grid_shape = [config.grid.nx, config.grid.ny, config.grid.nz]
span_average_shape = [config.grid.nx, config.grid.ny]

numwrites = int(config.temporal.num_iter / config.temporal.span_average_io_steps)
velocity_dset  = flowfields.create_dataset("velocity", 1, [5, *grid_shape], 1)
span_average_dset = flowfields.create_dataset("span_average", numwrites, [5, *span_average_shape], 1)

# now the main solver loop
streams.wrap_init_solver()

for i in range(config.temporal.num_iter):
    #print(f"step {i}")
    streams.wrap_step_solver()

    if (i % config.temporal.span_average_io_steps) == 0:
        hprint("writing span average to output")
        streams.wrap_copy_gpu_to_cpu()
        streams_data_slice = config.slice_flowfield_array(streams.mod_streams.w)
        calculate_span_averages(span_average, temp_field, streams_data_slice)

        span_average_dset.write_array(span_average)

# writes the full array out to as the final exit condition
streams.wrap_copy_gpu_to_cpu()
velocity_dset.write_array(config.slice_flowfield_array(streams.mod_streams.w))


#
# wrap up execution of solver
#

streams.wrap_finalize_solver()

print("finalizing solver")
streams.wrap_finalize()
