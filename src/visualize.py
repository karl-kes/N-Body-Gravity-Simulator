import struct, re, argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import OrderedDict
from pathlib import Path

# Usage: python visualize.py
#        python visualize.py --sim validation/sim_output.bin --speed 2

#   Space  = play / pause
#   Right  = step forward
#   Left   = step backward
#   Up     = speed up (2×)
#   Down   = slow down (0.5×)
#   R      = reset to frame 0
#   T      = toggle trails
#   Q      = quit


AU = 1.496e11

PLANETS = {"Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"}

COLORS = {
    "Sun": "#FDB813", "Mercury": "#8C7853", "Venus": "#E6C35C",
    "Earth": "#4A90D9", "Moon": "#AAAAAA", "Mars": "#C1440E",
    "Jupiter": "#C88B3A", "Saturn": "#E8D191", "Uranus": "#73C2C6",
    "Neptune": "#3B5BA5", "Pluto": "#D2B48C",
}

MARKER_SIZES = {
    "Sun": 120, "Jupiter": 70, "Saturn": 60, "Uranus": 45,
    "Neptune": 45, "Earth": 35, "Venus": 35, "Mars": 30,
    "Mercury": 20, "Pluto": 15,
}


def load(path):
    data = OrderedDict()
    with open(path, "rb") as f:
        # Read header
        num_bodies = struct.unpack("Q", f.read(8))[0]
        names = []
        for _ in range(num_bodies):
            raw = f.read(32)
            names.append(raw.split(b'\x00')[0].decode())

        # Initialize storage
        for name in names:
            data[name] = {"x": [], "y": [], "z": []}

        # Read frames
        doubles_per_frame = 2 + num_bodies * 6
        frame_size = doubles_per_frame * 8

        while True:
            raw = f.read(frame_size)
            if len(raw) < frame_size:
                break

            frame = np.frombuffer(raw, dtype=np.float64)
            states = frame[2:].reshape(num_bodies, 6)

            for i, name in enumerate(names):
                data[name]["x"].append(states[i, 0] / AU)
                data[name]["y"].append(states[i, 1] / AU)
                data[name]["z"].append(states[i, 2] / AU)

    # Convert lists to numpy arrays
    for name in data:
        data[name] = {k: np.array(v) for k, v in data[name].items()}

    return data


def read_config(path="src/Config.hpp"):
    try:
        text = Path(path).read_text()
        m = re.search(r'num_years\{\s*(\d+)\s*\}', text)
        return int(m.group(1)) if m else 100
    except FileNotFoundError:
        return 100


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sim", default="src/validation/sim_output.bin")
    p.add_argument("--speed", type=float, default=1.0)
    args = p.parse_args()

    data = load(args.sim)
    num_years = read_config()
    names = list(data.keys())
    n_frames = len(data[names[0]]["x"])

    print(f"Loaded {len(names)} bodies, {n_frames} epochs, {num_years} years")
    print("Controls: Space=play/pause, Arrows=step/speed, R=reset, T=trails, Q=quit")

    # Figure Setup:

    fig = plt.figure(figsize=(12, 10), facecolor="black")
    fig.canvas.manager.set_window_title("N-Body Solar System")
    ax = fig.add_subplot(111, projection="3d", facecolor="black")

    # Axis styling
    ax.set_xlabel("X (AU)", color="#555", fontsize=9, labelpad=8)
    ax.set_ylabel("Y (AU)", color="#555", fontsize=9, labelpad=8)
    ax.set_zlabel("Z (AU)", color="#555", fontsize=9, labelpad=8)
    ax.tick_params(colors="#444", labelsize=7)
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("#222")
    ax.yaxis.pane.set_edgecolor("#222")
    ax.zaxis.pane.set_edgecolor("#222")
    ax.grid(True, alpha=0.08, color="#666")

    # Compute axis limits from planet data only
    all_x, all_y, all_z = [], [], []
    for n in names:
        if n in PLANETS:
            all_x.extend(data[n]["x"])
            all_y.extend(data[n]["y"])
            all_z.extend(data[n]["z"])
    pad = 1.05
    lim = max(abs(min(all_x)), abs(max(all_x)),
              abs(min(all_y)), abs(max(all_y))) * pad
    zlim = max(abs(min(all_z)), abs(max(all_z))) * pad * 2
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-zlim, zlim)

    # Plot Elements:

    trails = {}
    dots = {}
    labels = {}

    for name in names:
        is_planet = name in PLANETS
        color = COLORS.get(name, "#666666")
        ms = MARKER_SIZES.get(name, 6 if not is_planet else 12)
        trail_alpha = 0.35 if is_planet else 0.1
        trail_lw = 1.0 if is_planet else 0.3

        trail, = ax.plot([], [], [], color=color, alpha=trail_alpha,
                         lw=trail_lw, zorder=1)
        dot = ax.scatter([], [], [], color=color, s=ms,
                         edgecolors="none", zorder=3, depthshade=True)
        trails[name] = trail
        dots[name] = dot

        if is_planet and name != "Sun":
            lbl = ax.text(0, 0, 0, f"  {name}", color="#777", fontsize=7,
                          fontfamily="monospace", ha="left", va="center", zorder=4)
            labels[name] = lbl

    title = ax.set_title("", color="#999", fontsize=11,
                         fontfamily="monospace", pad=12)

    # Animation:

    state = {"frame": 0, "playing": True, "speed": args.speed, "trails_on": True}

    def update(frame_num):
        f = state["frame"]

        for name in names:
            d = data[name]
            is_planet = name in PLANETS

            # Current position
            x, y, z = d["x"][f], d["y"][f], d["z"][f]
            dots[name]._offsets3d = ([x], [y], [z])

            # Trail
            if state["trails_on"]:
                trails[name].set_data_3d(d["x"][:f+1], d["y"][:f+1], d["z"][:f+1])
                trails[name].set_visible(True)
            else:
                trails[name].set_visible(False)

            # Label
            if name in labels:
                labels[name].set_position_3d((x, y, z))

        year = f / max(n_frames - 1, 1) * num_years
        title.set_text(f"Year {year:.1f} / {num_years}   |   "
                       f"{'▶' if state['playing'] else '⏸'}  {state['speed']}×   |   "
                       f"{len(names)} bodies")

        if state["playing"]:
            state["frame"] = min(f + 1, n_frames - 1)
            if state["frame"] >= n_frames - 1:
                state["frame"] = 0

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
        elif event.key == "t":
            state["trails_on"] = not state["trails_on"]
        elif event.key == "q":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)

    # Interval adapts to speed:
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