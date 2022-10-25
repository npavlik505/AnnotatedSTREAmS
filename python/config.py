import json
from typing import Any, Dict, Optional
import numpy as np

class Length():
    def __init__(self, lx: float, ly: float, lz: float):
        self.lx = lx
        self.ly = ly
        self.lz = lz

    @staticmethod
    def from_json(json_config: Dict[str, Any]):
        lx = json_config["x_length"]
        ly = json_config["y_length"]
        lz = json_config["z_length"]

        return Length(lx, ly, lz)

class Grid():
    def __init__(self, nx: int, ny: int, nz: int, ny_wr:int, ly_wr:float, dy_w:float, jb_grid:float, ng: int):
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.ny_wr = ny_wr
        self.ly_wr = ly_wr 

        self.dy_w = dy_w
        self.jb_gri = jb_grid

        self.ng = ng

    @staticmethod
    def from_json(json_config: Dict[str, Any]):
        nx = json_config["x_divisions"]
        ny = json_config["y_divisions"]
        nz = json_config["z_divisions"]

        ny_wr = 201
        ly_wr = 2.5
        dy_w = 0.7
        jb_grid = 0

        ng = 3

        return Grid(nx,ny,nz,ny_wr, ly_wr, dy_w, jb_grid, ng)

class Mpi():
    def __init__(self, x_split: int, z_split: int):
        self.x_split = x_split
        self.z_split = z_split

    @staticmethod
    def from_json(json_config: Dict[str, Any]):
        x_split = json_config["mpi_x_split"]
        z_split = 1

        return Mpi(x_split, z_split)

class Temporal():
    def __init__(self, num_iter: int, cfl: float, dt_control: int, print_control: int, io_type: int, span_average_io_steps: int):
        self.num_iter     = num_iter 
        self.cfl          = cfl
        self.dt_control   = dt_control 
        self.print_control= print_control
        self.io_type      = io_type

        # number of steps between writing writing the span average of all the rho / rhou / rhov / rhow/ energy
        self.span_average_io_steps = span_average_io_steps

    @staticmethod
    def from_json(json_config: Dict[str, Any]):
        num_iter = json_config["steps"]
        cfl = 0.75
        dt_control = 1
        print_control = 1
        io_type = 2

        span_average_io_steps = json_config["span_average_io_steps"]

        return Temporal(num_iter, cfl, dt_control, print_control, io_type, span_average_io_steps)

class Physics():
    def __init__(self, mach:float, reynolds_friction:float, temp_ratio: float, visc_type:int, Tref: float, turb_inflow: float ):
        self.mach             = mach
        self.reynolds_friction= reynolds_friction
        self.temp_ratio       = temp_ratio
        self.visc_type        = visc_type 
        self.Tref             = Tref
        self.turb_inflow      = turb_inflow

    @staticmethod
    def from_json(json_config: Dict[str, Any]):
        mach = json_config["mach_number"]
        reynolds_friction = json_config["reynolds_number"]
        temp_ratio = 1.
        visc_type = 2
        Tref = 160.
        turb_inflow = 0.75

        return Physics(mach, reynolds_friction, temp_ratio, visc_type, Tref, turb_inflow)

class Config():
    def __init__(self, length: Length, grid: Grid, mpi: Mpi, temporal: Temporal, physics: Physics):
        self.length = length
        self.grid = grid
        self.mpi = mpi
        self.temporal = temporal
        self.physics = physics

    @staticmethod
    def from_json(json_config: Dict[str, Any]):
        length = Length.from_json(json_config)
        grid = Grid.from_json(json_config)
        mpi = Mpi.from_json(json_config)
        temporal = Temporal.from_json(json_config)
        physics = Physics.from_json(json_config)

        return Config(length, grid, mpi, temporal, physics)

    def x_start(self) -> int:
        return self.grid.ng

    def x_end(self) -> int:
        return self.x_start() + int(self.grid.nx / self.mpi.x_split)

    def nx_mpi(self) -> int:
        return self.x_end() - self.x_start()

    def y_start(self) -> int:
        return self.grid.ng

    def y_end(self) -> int:
        return self.y_start() + self.grid.ny

    def ny_mpi(self) -> int:
        return self.y_end() - self.y_start()

    def z_start(self) -> int:
        return self.grid.ng

    def z_end(self) -> int:
        return self.z_start() + self.grid.nz

    def nz_mpi(self) -> int:
        return self.z_end() - self.z_start()

    def slice_flowfield_array(self, array: np.ndarray) -> np.ndarray:
        return array[:, \
                self.x_start():self.x_end(), \
                self.y_start():self.y_end(), \
                self.z_start():self.z_end()\
        ]
