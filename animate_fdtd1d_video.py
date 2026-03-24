"""Generate a video/GIF animation of 1D FDTD wave evolution.

This script runs a simple 1D FDTD solver (from fdtd1d.py) and saves an animation
showing the time evolution of the fields.

Usage:
    python animate_fdtd1d_video.py

Options:
    --output <file>        Output filename (ends with .gif or .mp4). Default: fdtd1d.gif
    --tfinal <float>       Simulation end time. Default: 0.5
    --frames <int>         Number of frames in the animation. Default: 120
    --boundary <b1> <b2>   Boundary condition strings, e.g. PEC PEC, PMC PMC, periodic periodic.
    --show                 Show animation window after saving (requires an interactive display).

The default run animates a Gaussian E-field pulse with PEC boundaries.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from fdtd1d import FDTD1D


def gaussian(x, x0, sigma):
    return np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def make_animation(
    x,
    initial_e,
    initial_h,
    boundaries,
    t_final=0.5,
    n_frames=120,
    outfile="fdtd1d.gif",
    show=False,
):
    """Create an animation of the E and H field evolution."""

    outfile = str((Path(__file__).resolve().parent / outfile).resolve())

    fdtd = FDTD1D(x, boundaries=boundaries)
    fdtd.load_initial_field(initial_e)
    if initial_h is not None:
        fdtd.load_initial_field(initial_h)

    # Precompute frame times and field snapshots
    times = np.linspace(0.0, t_final, n_frames)

    # Start from initial condition
    frames_e = [fdtd.get_e().copy()]
    frames_h = [fdtd.get_h().copy()]

    for t in times[1:]:
        fdtd.run_until(t)
        frames_e.append(fdtd.get_e().copy())
        frames_h.append(fdtd.get_h().copy())

    # Build animation
    fig, (ax_e, ax_h) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    line_e, = ax_e.plot(x, frames_e[0], lw=2, color="tab:blue")
    ax_e.set_ylabel("E(x)")
    ax_e.set_title("Electric field")
    ax_e.grid(alpha=0.3)

    line_h, = ax_h.plot(x[:-1], frames_h[0], lw=2, color="tab:orange")
    ax_h.set_ylabel("H(x)")
    ax_h.set_xlabel("x")
    ax_h.set_title("Magnetic field")
    ax_h.grid(alpha=0.3)

    # Set axis limits based on initial field magnitude
    e_max = np.max(np.abs(frames_e))
    h_max = np.max(np.abs(frames_h))
    ax_e.set_ylim(-1.1 * e_max, 1.1 * e_max)
    ax_h.set_ylim(-1.1 * h_max, 1.1 * h_max)

    time_txt = fig.suptitle("")

    def update(i):
        line_e.set_data(x, frames_e[i])
        line_h.set_data(x[:-1], frames_h[i])
        time_txt.set_text(f"t = {times[i]:.3f}")
        return line_e, line_h, time_txt

    anim = FuncAnimation(fig, update, frames=len(times), blit=True)

    # Choose writer based on file extension
    ext = Path(outfile).suffix.lower()
    if ext == ".gif":
        from matplotlib.animation import PillowWriter

        writer = PillowWriter(fps=25)
        anim.save(outfile, writer=writer)
    elif ext in (".mp4", ".m4v"):
        anim.save(outfile, fps=25, extra_args=["-vcodec", "libx264"])
    else:
        raise ValueError("Unsupported output extension: {}".format(ext))

    if show:
        plt.show()

    return outfile


def main():
    parser = argparse.ArgumentParser(description="Animate 1D FDTD evolution")
    parser.add_argument("--output", default="fdtd1d.gif", help="Output animation file")
    parser.add_argument("--tfinal", type=float, default=0.5, help="Final simulation time")
    parser.add_argument("--frames", type=int, default=120, help="Number of animation frames")
    parser.add_argument(
        "--boundary",
        nargs=2,
        default=["PEC", "PEC"],
        help="Boundary conditions (e.g. PEC PEC, PMC PMC, periodic periodic)",
    )
    parser.add_argument("--show", action="store_true", help="Show animation window after saving")
    args = parser.parse_args()

    x = np.linspace(-1.0, 1.0, 401)
    x0 = 0.0
    sigma = 0.08

    initial_e = gaussian(x, x0, sigma)
    initial_h = np.zeros_like(x[:-1])

    out = make_animation(
        x=x,
        initial_e=initial_e,
        initial_h=initial_h,
        boundaries=tuple(args.boundary),
        t_final=args.tfinal,
        n_frames=args.frames,
        outfile=args.output,
        show=args.show,
    )
    print(f"Saved animation: {out}")


if __name__ == "__main__":
    main()
