"""Create an animation of the 1D wave propagation solved by FDTD1D.

Run this script to generate an animation file (GIF by default) showing how the
field evolves in time.

Usage:
    python animate_fdtd1d.py

The script will create `fdtd1d.gif` in the current folder and optionally show
an interactive matplotlib window.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from fdtd1d import FDTD1D


from pathlib import Path


def make_animation(
    x,
    initial_field,
    t_final=0.5,
    n_frames=100,
    outfile="fdtd1d.gif",
    show=False,
):
    """Create and optionally display an animation of the wave evolution."""

    # Default to placing output next to this script to avoid permission issues.
    outfile = str((Path(__file__).resolve().parent / outfile).resolve())

    fdtd = FDTD1D(x)
    fdtd.load_initial_field(initial_field)

    fig, ax = plt.subplots()
    line, = ax.plot(x, initial_field, lw=2)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(initial_field.min() * 1.2, initial_field.max() * 1.2)
    ax.set_xlabel("x")
    ax.set_ylabel("E(x, t)")
    title = ax.set_title("")

    times = np.linspace(0.0, t_final, n_frames)

    def update(frame_index):
        t = times[frame_index]
        fdtd.run_until(t)
        e = fdtd.get_e()
        line.set_data(x, e)
        title.set_text(f"t = {t:.3f}")
        return line, title

    anim = FuncAnimation(fig, update, frames=len(times), blit=True)

    # Save as GIF; this requires Pillow (usually installed with matplotlib).
    writer = PillowWriter(fps=20)
    anim.save(outfile, writer=writer)

    if show:
        plt.show()

    return outfile


if __name__ == "__main__":
    # Example usage: Gaussian pulse
    x = np.linspace(-1.0, 1.0, 201)
    x0 = 0.0
    sigma = 0.05
    # Use a proper Gaussian peak (decays away from center)
    initial_e = np.exp(-0.5 * ((x - x0) / sigma) ** 2)

    out = make_animation(x, initial_e, t_final=0.3, n_frames=120, outfile="fdtd1d.gif", show=True)
    print(f"Saved animation: {out}")
