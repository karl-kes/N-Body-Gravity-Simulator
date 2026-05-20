"""
MW-Andromeda collision viewer (matplotlib).

Reads tests/sim_output.bin and displays the collision as a 3D scatter plot
with playback controls. First N/2 particles drawn blue (MW), last N/2 red (M31).

Usage:
    python src/visualize_galaxies.py
    python src/visualize_galaxies.py --speed 4

Controls:
    Space  play / pause
    Right  step forward
    Left   step backward
    Up     speed up (2x)
    Down   slow down (0.5x)
    R      reset to frame 0
    Q      quit
"""

import struct, argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

KPC = 3.086e19
GYR = 3.156e16


def load(path):
    """Load positions only from sim_output.bin.

    Returns (positions, times, N) where positions is (n_frames, N, 3) in kpc.
    """
    with open(path, "rb") as f:
        N = struct.unpack("Q", f.read(8))[0]
        f.read(N * 32)  # skip 32-byte name records

        doubles_per_frame = 2 + N * 6
        frame_size = doubles_per_frame * 8

        frames = []
        times = []
        while True:
            raw = f.read(frame_size)
            if len(raw) < frame_size:
                break
            arr = np.frombuffer(raw, dtype=np.float64)
            times.append(arr[1])
            states = arr[2:].reshape(N, 6)
            frames.append(states[:, 0:3] / KPC)

    return np.array(frames), np.array(times), N


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", default="tests/sim_output.bin")
    ap.add_argument("--speed", type=float, default=1.0)
    args = ap.parse_args()

    print(f"Loading {args.sim}...")
    positions, times, N = load(args.sim)
    n_frames = len(positions)
    half = N // 2

    print(f"  {N} particles, {n_frames} frames")
    print(f"  Time: {times[0] / GYR:.2f} Gyr -> {times[-1] / GYR:.2f} Gyr")
    print("Controls: Space=play/pause, Arrows=step/speed, R=reset, Q=quit")

    fig = plt.figure(figsize=(12, 10), facecolor="black")
    fig.canvas.manager.set_window_title("MW-Andromeda Collision")
    ax = fig.add_subplot(111, projection="3d", facecolor="black")

    ax.set_xlabel("X (kpc)", color="#555", fontsize=9, labelpad=8)
    ax.set_ylabel("Y (kpc)", color="#555", fontsize=9, labelpad=8)
    ax.set_zlabel("Z (kpc)", color="#555", fontsize=9, labelpad=8)
    ax.tick_params(colors="#444", labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor("#222")
    ax.grid(True, alpha=0.08, color="#666")

    # Axis limits from the full trajectory (with 5% padding).
    lim = np.abs(positions).max() * 1.05
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    # Two scatter objects: MW (blue) and M31 (red).
    # Small markers, low alpha so dense regions show structure without saturating.
    mw = ax.scatter(
        positions[0, :half, 0], positions[0, :half, 1], positions[0, :half, 2],
        c="#4A90D9", s=0.4, alpha=0.6, edgecolors="none", depthshade=False,
    )
    m31 = ax.scatter(
        positions[0, half:, 0], positions[0, half:, 1], positions[0, half:, 2],
        c="#E64545", s=0.4, alpha=0.6, edgecolors="none", depthshade=False,
    )
    title = ax.set_title("", color="#999", fontsize=11, fontfamily="monospace", pad=12)

    state = {"frame": 0, "playing": True, "speed": args.speed}

    def update(_):
        f = state["frame"]
        mw._offsets3d  = (positions[f, :half, 0],
                          positions[f, :half, 1],
                          positions[f, :half, 2])
        m31._offsets3d = (positions[f, half:, 0],
                          positions[f, half:, 1],
                          positions[f, half:, 2])
        t_gyr = times[f] / GYR
        title.set_text(f"t = {t_gyr:.2f} Gyr   |   "
                       f"{'>' if state['playing'] else '||'}  {state['speed']}x   |   "
                       f"{N} particles  (MW + M31)")
        return []

    def on_key(event):
        if event.key == " ":
            state["playing"] = not state["playing"]
        elif event.key == "right":
            state["frame"] = min(state["frame"] + 1, n_frames - 1)
        elif event.key == "left":
            state["frame"] = max(state["frame"] - 1, 0)
        elif event.key == "up":
            state["speed"] = min(state["speed"] * 2, 32)
        elif event.key == "down":
            state["speed"] = max(state["speed"] / 2, 0.25)
        elif event.key == "r":
            state["frame"] = 0
        elif event.key == "q":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)

    def frame_gen():
        while True:
            if state["playing"]:
                step = max(1, int(state["speed"]))
                state["frame"] = (state["frame"] + step) % n_frames
            yield state["frame"]

    anim = FuncAnimation(fig, update, frames=frame_gen, interval=33,
                         blit=False, cache_frame_data=False)

    plt.subplots_adjust(left=0, right=1, top=0.95, bottom=0.02)
    plt.show()


if __name__ == "__main__":
    main()
