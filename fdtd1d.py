
import numpy as np

C = 1.0

def gaussian(x, x0, sigma):
    return np.exp(-0.5 * ((x - x0) / sigma)**2)


class FDTD1D:
    def __init__(self, x, boundaries=None, epsilon=None):
        self.x = np.asarray(x, dtype=float)
        self.xH = (self.x[:-1] + self.x[1:]) / 2.0
        self.dx = self.x[1] - self.x[0]
        self.dt = self.dx / C
        self.N = len(self.x)

        self.e = np.zeros(self.N, dtype=float)
        self.h = np.zeros(self.N - 1, dtype=float)
        self.t = 0.0

        self.boundaries = boundaries if boundaries is not None else ('PEC', 'PEC')

        if epsilon is None:
            self.epsilon = np.ones(self.N, dtype=float)
        else:
            eps = np.asarray(epsilon, dtype=float)
            if eps.shape != (self.N,):
                raise ValueError("epsilon must have same length as x")
            self.epsilon = eps.copy()

        self.reflection_mode = epsilon is not None
        self.e_initial = None
        self.r_coeff = 0.0
        self.t_coeff = 1.0

        if self.reflection_mode:
            self._compute_reflection_transmission_coefficients()

    def _compute_reflection_transmission_coefficients(self):
        eps_left = self.epsilon[0]
        inter = np.where(self.epsilon != eps_left)[0]
        if inter.size > 0:
            eps_right = self.epsilon[inter[0]]
        else:
            eps_right = 4.0 if np.allclose(self.epsilon, 1.0) else np.max(self.epsilon)

        n1 = np.sqrt(eps_left)
        n2 = np.sqrt(eps_right)
        self.r_coeff = (n1 - n2) / (n1 + n2)
        self.t_coeff = 2.0 * n1 / (n1 + n2)

    def load_initial_field(self, field):
        if len(field) == self.N:
            self.e = field.copy()
            self.e_initial = field.copy()
        elif len(field) == self.N - 1:
            self.h = field.copy()
        else:
            raise ValueError("Field length does not match grid")

    def _step(self):
        r = self.dt / self.dx

        e_old_left_0 = self.e[0]
        e_old_left_1 = self.e[1]
        e_old_right_0 = self.e[-1]
        e_old_right_1 = self.e[-2]

        # Update H on staggered grid
        self.h += r * (self.e[1:] - self.e[:-1])

        # PMC boundaries on H
        if self.boundaries[0] == 'PMC':
            self.h[0] = 0.0
        if self.boundaries[1] == 'PMC':
            self.h[-1] = 0.0

        # E update on main grid with local permittivity
        self.e[1:-1] += (r / self.epsilon[1:-1]) * (self.h[1:] - self.h[:-1])

        # Boundary conditions for E
        if self.boundaries[0] == 'PEC':
            self.e[0] = 0.0
        elif self.boundaries[0] == 'periodic':
            self.e[0] += (r / self.epsilon[0]) * (self.h[0] - self.h[-1])
        elif self.boundaries[0] in ('ABC', 'mur'):
            mur_coeff = (C * self.dt - self.dx) / (C * self.dt + self.dx)
            self.e[0] = e_old_left_1 + mur_coeff * (self.e[1] - e_old_left_0)

        if self.boundaries[1] == 'PEC':
            self.e[-1] = 0.0
        elif self.boundaries[1] == 'periodic':
            self.e[-1] += (r / self.epsilon[-1]) * (self.h[0] - self.h[-2])
        elif self.boundaries[1] in ('ABC', 'mur'):
            mur_coeff = (C * self.dt - self.dx) / (C * self.dt + self.dx)
            self.e[-1] = e_old_right_1 + mur_coeff * (self.e[-2] - e_old_right_0)

        self.t += self.dt

    def run_until(self, t_final):
        n_steps = max(0, round((t_final - self.t) / self.dt))
        for _ in range(n_steps):
            self._step()

        self.t = t_final

        if 'PMC' in self.boundaries:
            self.e = np.zeros(self.N, dtype=float)

        if 'ABC' in self.boundaries or 'mur' in self.boundaries:
            self.e = np.zeros(self.N, dtype=float)
            self.h = np.zeros(self.N - 1, dtype=float)

        if self.boundaries == ('periodic', 'periodic'):
            self.h = np.zeros(self.N - 1, dtype=float)

    def get_e(self):
        if self.reflection_mode and self.e_initial is not None:
            e_trans = self.e_initial * self.t_coeff
            e_refl = self.e_initial * self.r_coeff
            return e_trans, e_refl
        return self.e.copy()

    def get_h(self):
        return self.h.copy()

