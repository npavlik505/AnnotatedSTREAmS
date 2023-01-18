import numpy as np
from config import Config
import libstreams as streams

# the equation of the polynomial for the jet actuation in coordinates local to the jet
# actuator. This means that the jet actuator starts at x = 0 and ends at x = slot_end
# 
# this must be recomputed at every amplitude change
class Polynomial():
    def __init__(self, a: float, b: float, c: float):
        self.a = a
        self.b = b
        self.c = c

    def evaluate(self, x_idx: int) -> float:
        return self.a * x_idx**2  + \
                self.b * x_idx + \
                self.c

# helper class to recompute the polynomial of the jet actuator
class PolynomialFactory():
    def __init__(self, vertex_x: float, slot_start: float, slot_end: float):
        self.slot_start = slot_start
        self.slot_end = slot_end
        self.vertex_x  = vertex_x

    def poly(self, amplitude: float) -> Polynomial:
        a = amplitude/(self.vertex_x**2 - self.vertex_x*self.slot_end- (self.vertex_x - self.slot_end)*self.slot_start)
        b = -(amplitude*self.slot_end+ amplitude*self.slot_start)/(self.vertex_x**2 - self.vertex_x*self.slot_end- (self.vertex_x - self.slot_end)*self.slot_start)
        c = amplitude*self.slot_end*self.slot_start/(self.vertex_x**2 - self.vertex_x*self.slot_end- (self.vertex_x - self.slot_end)*self.slot_start)

        return Polynomial(a, b, c)

class JetActuator():
    def __init__(self, rank: int, config: Config):
        # if x_start_slot == -1 then we do not have a matrix
        # allocated on this MPI process

        self.rank = rank
        self.config = config

        vertex_x = (config.jet.slot_start +  config.jet.slot_end) / 2
        self.factory = PolynomialFactory(vertex_x, config.jet.slot_start, config.jet.slot_end)

        self.local_slot_start_x = streams.mod_streams.x_start_slot
        self.local_slot_nx = streams.mod_streams.nx_slot
        self.local_slot_nz = streams.mod_streams.nz_slot

        self.has_slot = streams.mod_streams.x_start_slot != -1

        if self.has_slot:
            #print(self.local_slot_nx, self.local_slot_nz, self.local_slot_start_x)
            self.bc_velocity = streams.mod_streams.blowing_bc_slot_velocity[:, :]

    def set_amplitude(self, amplitude: float):
        # WARNING: copying to GPU and copying to CPU must happen on ALL mpi procs
        # if we have a matrix, we can adjust the velocity of the blowing actuator
        streams.wrap_copy_blowing_bc_to_cpu()

        if not self.has_slot:
            streams.wrap_copy_blowing_bc_to_gpu()
            return None

        # calculate the equation of the polynomial that the velocity of the jet
        # actuator will have with this amplitude
        poly = self.factory.poly(amplitude)

        for idx in range(self.local_slot_nx):
            local_x = self.local_slot_start_x + idx
            global_x = self.config.local_to_global_x(local_x, self.rank)

            velo =  poly.evaluate(global_x)
            #print(f"velocity at idx {idx} / local x {local_x} / global x {global_x} ::: {velo}")

            self.bc_velocity[idx, 0:self.local_slot_nz] = velo

        #print(f"array shape for BC", self.bc_velocity.shape)

        # finally, copy everything back to the GPU
        streams.wrap_copy_blowing_bc_to_gpu()
        return None

