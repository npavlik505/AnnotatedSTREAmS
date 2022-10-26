#
# Mostly harvested from
# https://github.com/Fluid-Dynamics-Group/selective-modification/blob/763f1b6369a851b374bb918270ac8a80f72f5738/solver/io_utils.py
#
import h5py
from mpi4py import MPI
import numpy as np
from typing import Tuple, Optional, List, Any

class Group():
    def __init__(self, group: h5py.Group, rank: int):
        self.group = group
        self.rank = rank

    def write_attr(self, attr: str, value: Any):
        self.group[attr] = value

class Dataset():
    def __init__(self, dset: h5py.Dataset, num_writes:int, array_shape: List[int], mpi_split_idx: int, rank: int):
        self.num_writes = num_writes
        self.step_number = 0
        self.split_idx = mpi_split_idx
        self.dim = len(array_shape)
        self.rank = rank

        self.dset = dset

    def write_array(self, array: np.ndarray):
        # ensure we are not overwriting on the indicies.
        # if we are, we will get a nicer error than h5 will provide
        if self.step_number >= self.num_writes:
            import warnings
            if MPI.COMM_WORLD.Get_rank() == 0:
                warnings.warn(f"attempted to write {self.step_number+1} values to velocity h5 file, when constructed for {self.num_writes} writes - skipping this write", RuntimeWarning)
            return None
        
        if self.dim != len(array.shape):
            raise ValueError(f"reported dimension ({self.dim}) in initialization is different from argument `array` {len(array.shape)}")

        # writing 3D vector fields (partial writes by each mpi process)
        if self.dim == 4 and self.split_idx == 1:
            # find out the position that we are writing to
            split_size = np.size(array,self.split_idx)
            start_slice = split_size * self.rank
            end_slice = start_slice + split_size

            # (write step number, vector components ex: (u,v,w,), x data, y data, z data)
            self.dset[self.step_number, :, start_slice:end_slice, :, :] = array

        # writing 3D scalar fields (partial writes by each mpi process)
        elif self.dim == 3 and self.split_idx == 1:
            # find out the position that we are writing to
            split_size = np.size(array,self.split_idx)
            start_slice = split_size * self.rank
            end_slice = start_slice + split_size

            #print(f"writing to {start_slice}:{end_slice}")

            self.dset[self.step_number, :, start_slice:end_slice, :] = array
        
        # writing a 1D array long the x axis
        elif self.dim == 1 and self.split_idx == 0:
            split_size = np.size(array,self.split_idx)
            start_slice = split_size * self.rank
            end_slice = start_slice + split_size

            self.dset[self.step_number, start_slice:end_slice] = array

        # writing spectra arrays
        elif self.dim == 2 and self.split_idx < 0:
            self.dset[self.step_number, :] = array

        # writing timesteps
        elif self.dim == 1 and self.split_idx < 0:
            self.dset[self.step_number, :] = array

        else:
            raise ValueError(f"unhandled dimension {self.dim} and mpi split along axis {self.split_idx}")

        self.step_number += 1;

# Write a h5 file using MPI
#
# for a numpy array of shape = (3, nx, ny, nz), with num_proc number of mpi processesses,
# this class will write a giant array to disk in a 6 dimensional array in the shape of 
# (num_writes, 3, nx, ny, nz)
class IoFile:
    # params
    # filename: the path to the .h5 file that you wish to write to
    # num_writes: the total number of times that you will call `write_array`
    def __init__(self, filename: str):
        comm = MPI.COMM_WORLD
        self.rank = comm.rank

        self.file = h5py.File(filename, 'w', driver='mpio', comm = MPI.COMM_WORLD)

    # mpi_split_idx: in the array you send to write to the HDF5 file, which 0-based index is the MPI data 
    #   divided upon? For example, streams is (ideally) portioned by MPI along the x-axis, so if you have a vector
    #   field (3, nx, ny, nz), then mpi_split_idx = 1 so that we tell the library that nx is where things are split upon
    def create_dataset(self, name: str, num_writes: int, array_shape: List[int], mpi_split_idx: int) -> Dataset:

        dataset = self.file.create_dataset(name, (num_writes, *array_shape))

        return Dataset(dataset, num_writes, array_shape, mpi_split_idx, self.rank)

    def create_group(self, name: str) -> Group:
        group = self.file.create_group(name)
        return Group(group, self.rank)

    # close the underlying h5 file handle
    def close(self):
        self.file.close()

def write_initial_condition(base_path: str, velocity: np.ndarray, packed_spec: np.ndarray):
    VALIDATE_IC_WRITING = False

    # statement for validating the writing of the initial condition
    if VALIDATE_IC_WRITING:
        array = np.zeros(velocity.shape)
        proc = MPI.COMM_WORLD.rank
        base_number = proc * 3
        array[0, :, :, :] = base_number
        array[1, :, :, :] = base_number + 1
        array[2, :, :, :] = base_number + 2

    vec, _nx, ny, nz = velocity.shape
    # we know what the shape should be here globally,
    # which is not reflected by array.shape

    shape = [vec, ny, ny, nz]

    h5_path = f"{base_path}/initial_condition.h5"
    f = IoFile(h5_path)

    velocity_dset = f.create_dataset("velocity", 1, shape, 1)
    velocity_dset.write_array(velocity)

    assert(len(packed_spec.shape) == 2)

    spectra_dset = f.create_dataset("spectra", 1, list(packed_spec.shape), -1)
    spectra_dset.write_array(packed_spec)

    f.close()

