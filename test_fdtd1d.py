import numpy as np
import matplotlib.pyplot as plt
import pytest
from fdtd1d import FDTD1D, C, gaussian




def test_fdtd_solves_basic_propagation():
    x = np.linspace(-1, 1, 201)
    x0 = 0.0
    sigma = 0.05
    initial_e = gaussian(x, x0, sigma)

    fdtd = FDTD1D(x)
    fdtd.load_initial_field(initial_e)

    t_final = 0.2
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()

    e_expected = 0.5 * gaussian(x, -t_final*C, sigma) \
     + 0.5 * gaussian(x, t_final*C, sigma)

    plt.plot(x, e_solved)
    plt.plot(x, e_expected)

    assert np.corrcoef(e_solved, e_expected)[0,1] > 0.99


def test_fdtd_PEC_boundary_conditions():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    boundaries = ('PEC', 'PEC')

    x0 = 0.0
    sigma = 0.05
    initial_e = gaussian(x, x0, sigma)

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)

    L = xMax - xMin
    t_final = L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()

    e_expected = - initial_e
    h_expected = np.zeros_like(h_solved)

    assert np.corrcoef(e_solved, e_expected)[0,1] > 0.99
    assert np.allclose(h_solved, h_expected, atol=0.01)

def test_fdtd_PMC_boundary_conditions():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    boundaries = ('PMC', 'PMC')
    
    x0 = 0.0
    sigma = 0.05
    # H is stored on staggered cells of size N-1
    x_h = x[:-1]
    initial_h = gaussian(x_h, x0, sigma)

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_h)

    L = xMax - xMin
    t_final = 0.5*L / C
    fdtd.run_until(t_final)

    h_solved = fdtd.get_h()
    e_solved = fdtd.get_e()

    h_expected = - initial_h
    e_expected = np.zeros_like(e_solved)
    
    assert np.corrcoef(h_solved, h_expected)[0,1] > 0.98
    assert np.allclose(e_solved, e_expected, atol=0.01)
    

def test_fdtd_periodic_boundary_conditions():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    boundaries = ('periodic', 'periodic')

    x0 = 0.0
    sigma = 0.05
    initial_e = gaussian(x, x0, sigma)
    initial_h = np.zeros_like(initial_e[:-1])

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)

    L = xMax - xMin
    t_final = L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()

    e_expected = initial_e
    h_expected = initial_h

    assert np.corrcoef(e_solved, e_expected)[0, 1] > 0.99
    assert np.allclose(h_solved, h_expected, atol=1e-10)

<<<<<<< HEAD
def test_fdtd_ABC_boundary_conditions():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    boundaries = ('PEC', 'ABC')
    
    x0 = 0.0
    sigma = 0.05
    initial_e = gaussian(x, x0, sigma)

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)

    L = xMax - xMin
    t_final = L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()

    e_expected = np.zeros_like(e_solved)
    h_expected = np.zeros_like(h_solved)
    
    assert np.allclose(e_solved, e_expected, atol=0.01)
    assert np.allclose(h_solved, h_expected, atol=0.01)
=======
def test_fdtd_PMC_boundary_conditions():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    xH = (x[1:] + x[:-1]) / 2.0
    boundaries = ('PMC', 'PMC')
    
    x0 = 0.0
    sigma = 0.05
    initial_h = gaussian(xH, x0, sigma)
    initial_e = np.zeros_like(x)
>>>>>>> upstream/main

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)
    fdtd.h = initial_h.copy()

    L = xMax - xMin
    t_final = L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()

    h_expected = - initial_h
    e_expected = np.zeros_like(e_solved)

    assert np.corrcoef(h_solved, h_expected)[0,1] > 0.99
    assert np.max(np.abs(e_solved - e_expected)) < 1e-8

def test_fdtd_total_spread_field():
    xMax = 1
    xMin = -1
    L = xMax - xMin

    x = np.linspace(xMin,xMax,201)
    xH = (x[1:] + x[:-1]) / 2.0
    boundaries = ('mur','mur') 
    x_o = 0
    initial_e = np.zeros_like(x)

    def my_pert(t):
        return np.sin(t*np.pi)

    fdtd = FDTD1D(x,boundaries,x_o,pert = lambda t: my_pert(t))
    fdtd.load_initial_field(initial_e)
    t_final = L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()

    e_expected = my_pert(x - L)
    h_expected = -my_pert(xH - L)
    
    assert np.corrcoef(e_solved, e_expected)[0,1] > -0.6
    assert np.corrcoef(h_solved, h_expected)[0,1] > 0.8

def test_fdtd_mur_boundary_conditions():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    xH = (x[1:] + x[:-1]) / 2.0
    boundaries = ('mur', 'mur')

    x0 = 0.0
    sigma = 0.05

    initial_e = gaussian(x, x0, sigma)
    initial_h = -gaussian(xH, x0, sigma)

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)
    fdtd.h = initial_h.copy()

    L = xMax - xMin
    t_final = 0.5*L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()

    assert np.allclose(e_solved, 0.0, atol=1e-2)
    assert np.allclose(h_solved, 0.0, atol=1e-2)
    
def test_reflection():
    xMax = 1
    xMin = -1
    x = np.linspace(xMin, xMax, 201)
    
    x0 = 0.0
    sigma = 0.05
    initial_e = gaussian(x, x0, sigma)
    
    L = xMax - xMin
    t_final = 2.0*L/C
    
    epsilon_r = 4.0
    
    epsilon = np.ones_like(x)
    epsilon[x > L] = epsilon_r
    
    n1 = 1.0
    n2 = np.sqrt(epsilon_r)
    
    R = (n1-n2)/(n1+n2)
    T = 2.0*n1/(n1+n2)
    
    fdtd = FDTD1D(x, epsilon=epsilon)
    fdtd.load_initial_field(initial_e)
    fdtd.run_until(t_final)
    
    e_trans_solved, e_reflex_solved = fdtd.get_e()
    
    e_reflex_expected = initial_e*R
    e_trans_expected = initial_e*T
    
    assert np.corrcoef(e_reflex_solved, e_reflex_expected)[0,1] > 0.99
    assert np.corrcoef(e_trans_solved, e_trans_expected)[0,1] > 0.99


from fdtd1d_dispersive import simulate_dispersive_slab_reflection_transmission


def test_dispersive_slab_reflection_transmission():
    x = np.linspace(-1.0, 1.0, 801)
    slab_start = -0.2
    slab_end = 0.2
    eps_inf = 1.0
    eps_s = 4.0
    tau = 0.15
    f0 = 2.0
    t_final = 6.0

    R_sim, T_sim, R_ana, T_ana = simulate_dispersive_slab_reflection_transmission(
        x,
        slab_start,
        slab_end,
        eps_inf,
        eps_s,
        tau,
        f0,
        t_final,
    )

    assert abs(R_sim - abs(R_ana)) < 0.20
    assert abs(T_sim - abs(T_ana)) < 0.20

def test_fdtd_dissipative_exact():
    xMax = 1.0
    xMin = -1.0
    L = xMax - xMin

    x = np.linspace(xMin, xMax, 201)
    xH = (x[1:] + x[:-1]) / 2.0
    boundaries = ('periodic', 'periodic')

    x0 = 0.0
    k = np.pi  
    initial_e = np.sin((x-x0)*k)
    initial_h = np.zeros_like(initial_e[:-1])

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)
    fdtd.h = initial_h.copy()
    
    sigma_val = 1.0
    eps_r_val = 2.0
    
    fdtd.sig = np.full_like(x, sigma_val)       
    fdtd.eps_r = np.full_like(x, eps_r_val)
    fdtd.eps = fdtd.eps0 * fdtd.eps_r

    gamma = sigma_val / (2 * fdtd.eps0 * eps_r_val)
    w0_sq = (k**2) / (fdtd.mu0 * fdtd.eps0 * eps_r_val)
    wd = np.sqrt(w0_sq - gamma**2)

    t_final = L / C
    fdtd.run_until(t_final)

    e_solved = fdtd.get_e()
    h_solved = fdtd.get_h()
    
    t_func_E = np.exp(-gamma * t_final) * (np.cos(wd * t_final) - (gamma / wd) * np.sin(wd * t_final))
    t_func_H = np.exp(-gamma * t_final) * (k / (fdtd.mu0 * wd)) * np.sin(wd * t_final)

    e_expected = t_func_E * np.sin(k * (x - x0))
    h_expected = -t_func_H * np.cos(k * (xH - x0)) 

    assert np.allclose(e_solved, e_expected, atol=1e-2)
    assert np.allclose(h_solved, h_expected, atol=1e-2)


def test_fdtd_dielectric_reflection():
    L = 2.0
    N = 401
    x = np.linspace(0, L, N)
    xH = (x[1:] + x[:-1]) / 2.0
    dx = x[1] - x[0]
    dt = dx / C
    
    boundaries = ('mur', 'mur')

    x0 = 0.4
    sigma = 0.05
    initial_e = gaussian(x, x0, sigma)
    initial_h = -gaussian(xH, x0, sigma)

    E_inc_max = np.max(initial_e)

    fdtd = FDTD1D(x, boundaries)
    fdtd.load_initial_field(initial_e)
    fdtd.h = initial_h.copy()

    eps_r_1 = 1.0   
    eps_r_2 = 4.0   
    interface_pos = L / 2
    
    fdtd.eps_r = np.where(x < interface_pos, eps_r_1, eps_r_2)
    fdtd.sig = np.zeros_like(x) 

    eta_0 = 1.0 / np.sqrt(eps_r_1)  # = 1
    eta_1 = 1.0 / np.sqrt(eps_r_2)  # = 0.5
    R_theory = (eta_0 - eta_1) / (eta_0 + eta_1)

    obs_idx = 30  # x ≈ 0.15

    t_at_interface = (interface_pos - x0) / C
    t_ref_at_obs = t_at_interface + (interface_pos - x[obs_idx]) / C

    t_final = 2.0
    n_steps = int(t_final / dt)
    
    E_ref_max_observed = 0.0
    
    for _ in range(n_steps):
        fdtd._step()
        if fdtd.t > t_ref_at_obs - 0.1:
            E_at_obs = np.abs(fdtd.e[obs_idx])
            E_ref_max_observed = max(E_ref_max_observed, E_at_obs)

    R_numerical = E_ref_max_observed / E_inc_max

    assert np.abs(R_numerical - np.abs(R_theory)) < 0.02


if __name__ == "__main__":
    pytest.main([__file__])