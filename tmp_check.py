import numpy as np
from fdtd1d import FDTD1D, C

x = np.linspace(-1, 1, 201)

# ABC case
fdtd = FDTD1D(x, boundaries=('PEC', 'ABC'))
initial_e = np.exp(-0.5 * ((x - 0.0) / 0.05) ** 2)
fdtd.load_initial_field(initial_e)
fdtd.run_until(2.0 / C)

es = fdtd.get_e()
hs = fdtd.get_h()
print('abc e max', np.max(np.abs(es)), 'h max', np.max(np.abs(hs)))
print('abc e nonzero', np.count_nonzero(np.abs(es) > 0.01), 'h nonzero', np.count_nonzero(np.abs(hs) > 0.01))

# PMC case
x_h = x[:-1]
fdtd = FDTD1D(x, boundaries=('PMC', 'PMC'))
initial_h = np.exp(-0.5 * ((x_h - 0.0) / 0.05) ** 2)
fdtd.load_initial_field(initial_h)
fdtd.run_until(2.0 / C)\n
hs2 = fdtd.get_h()
es2 = fdtd.get_e()
print('pmc h corr', np.corrcoef(hs2, -initial_h)[0, 1], 'e max', np.max(np.abs(es2)))
