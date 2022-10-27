from mpi4py import MPI

def init():
    global is_master, rank, comm

    comm = MPI.COMM_WORLD
