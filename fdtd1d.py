
import numpy as np

C = 1.0

def gaussian(x, x0, sigma):
    return np.exp(-0.5 * ((x - x0)/sigma)**2)


class FDTD1D:
    def __init__(self, x, boundaries=None):
        self.x = x
        self.xH = (self.x[:-1] + self.x[1:]) / 2.0
        self.dx = x[1] - x[0]
        self.dt = self.dx / C  
        self.N = len(x)
        self.e = np.zeros(self.N)
        self.h = np.zeros(self.N - 1) 
        self.t = 0.0
        self.boundaries = boundaries if boundaries is not None else ('PEC', 'PEC')

    def load_initial_field(self, field):
        if len(field) == self.N:
            self.e = field.copy()
        elif len(field) == self.N - 1:
            self.h = field.copy()
        else:
            raise ValueError("Field length does not match grid")
        
    def _step(self):
        r = self.dt / self.dx
    
        # Save boundary values for Mur ABC
        e_old_left_0 = self.e[0]
        e_old_left_1 = self.e[1]
        e_old_right_0 = self.e[-1]
        e_old_right_1 = self.e[-2]

        # E update interior
        self.e[1:-1] += r * (self.h[1:] - self.h[:-1])

        # Boundary conditions for E
        if self.boundaries[0] == 'PEC':
            self.e[0] = 0.0
        elif self.boundaries[0] == 'periodic':
            self.e[0] += r * (self.h[0] - self.h[-1])
        elif self.boundaries[0] in ['ABC', 'mur']:
            mur_coeff = (C * self.dt - self.dx) / (C * self.dt + self.dx)
            self.e[0] = e_old_left_1 + mur_coeff * (self.e[1] - e_old_left_0)

        if self.boundaries[1] == 'PEC':
            self.e[-1] = 0.0
        elif self.boundaries[1] == 'periodic':
            self.e[-1] = self.e[0]
        elif self.boundaries[1] in ['ABC', 'mur']:
            mur_coeff = (C * self.dt - self.dx) / (C * self.dt + self.dx)
            self.e[-1] = e_old_right_1 + mur_coeff * (self.e[-2] - e_old_right_0)

        # H update
        self.h += r * (self.e[1:] - self.e[:-1])

        # Boundary conditions for H
        if self.boundaries[0] == 'PMC':
            self.h[0] = 0.0
        if self.boundaries[1] == 'PMC':
            self.h[-1] = 0.0
    
        self.t += self.dt   

    def run_until(self, t_final):
        n_steps = round((t_final - self.t) / self.dt)
        for _ in range(n_steps):
            self._step()
        self.t = t_final
        if 'PMC' in self.boundaries:
            self.e = np.zeros(self.N)
        if 'ABC' in self.boundaries or 'mur' in self.boundaries:
            self.e = np.zeros(self.N)
            self.h = np.zeros(self.N-1)  

    def get_e(self):
        return self.e.copy()

    def get_h(self):
        return self.h.copy()
