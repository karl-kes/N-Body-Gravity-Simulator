import argparse, csv, json, re, struct, sys, time, urllib.request, urllib.parse
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

# Usage (reads config from Config.hpp automatically):
# python jpl_compare.py fetch --moons
# python jpl_compare.py compare

AU_KM = 1.496e8
SEC_PER_YR = 365.25 * 86400.0
G_KM3 = 6.6743e-20

def read_config(path="src/Config.hpp"):
    text = Path(path).read_text()
    cfg = {}
    for name in ["output_hours", "num_years"]:
        m = re.search(rf'{name}\{{\s*(\d+)\s*\}}', text)
        if not m:
            print(f"ERROR: Could not find '{name}' in {path}")
            sys.exit(1)
        cfg[name] = int(m.group(1))
    return cfg

PLANETS = [
    ("10",  "Sun",     1.32712440041e11, None),
    ("199", "Mercury", 2.2032e4,         None),
    ("299", "Venus",   3.24859e5,        None),
    ("399", "Earth",   3.98600e5,        None),
    ("499", "Mars",    4.28284e4,        None),
    ("599", "Jupiter", 1.26687e8,        None),
    ("699", "Saturn",  3.79312e7,        None),
    ("799", "Uranus",  5.79394e6,        None),
    ("899", "Neptune", 6.83653e6,        None),
    ("999", "Pluto",   8.71e2,           None),
]

MOONS = [
    ("301", "Moon",      4.9028e3,   "Earth"),
    ("401", "Phobos",    7.093e-4,   "Mars"),
    ("402", "Deimos",    9.615e-5,   "Mars"),
    ("501", "Io",        5.9600e3,   "Jupiter"),
    ("502", "Europa",    3.2027e3,   "Jupiter"),
    ("503", "Ganymede",  9.8878e3,   "Jupiter"),
    ("504", "Callisto",  7.1793e3,   "Jupiter"),
    ("505", "Amalthea",  1.38e-1,    "Jupiter"),
    ("601", "Mimas",     2.503e0,    "Saturn"),
    ("602", "Enceladus", 7.211e0,    "Saturn"),
    ("603", "Tethys",    4.121e1,    "Saturn"),
    ("604", "Dione",     7.312e1,    "Saturn"),
    ("605", "Rhea",      1.539e2,    "Saturn"),
    ("606", "Titan",     8.978e3,    "Saturn"),
    ("607", "Hyperion",  3.718e-1,   "Saturn"),
    ("608", "Iapetus",   1.205e2,    "Saturn"),
    ("609", "Phoebe",    5.531e-1,   "Saturn"),
    ("701", "Ariel",     8.346e1,    "Uranus"),
    ("702", "Umbriel",   8.509e1,    "Uranus"),
    ("703", "Titania",   2.269e2,    "Uranus"),
    ("704", "Oberon",    2.053e2,    "Uranus"),
    ("705", "Miranda",   4.4e0,      "Uranus"),
    ("801", "Triton",    1.4276e3,   "Neptune"),
    ("802", "Nereid",    2.06e0,     "Neptune"),
    ("901", "Charon",    1.0588e2,   "Pluto"),
]

# Horizons API:

API = "https://ssd.jpl.nasa.gov/api/horizons.api"

def fetch_vectors(body_id, start, stop, step):
    params = {
        "format": "text", "COMMAND": f"'{body_id}'",
        "OBJ_DATA": "NO", "MAKE_EPHEM": "YES", "EPHEM_TYPE": "VECTORS",
        "CENTER": "'500@0'", "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'",
        "STEP_SIZE": f"'{step}'", "REF_PLANE": "ECLIPTIC", "REF_SYSTEM": "ICRF",
        "VEC_TABLE": "2", "OUT_UNITS": "KM-S", "VEC_LABELS": "YES", "CSV_FORMAT": "YES",
    }
    url = API + "?" + urllib.parse.urlencode(params, safe="'@")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception as e:
            if attempt < 2: time.sleep(2*(attempt+1))
            else: raise RuntimeError(f"Failed {body_id}: {e}")

def parse_vectors(text):
    soe, eoe = text.find("$$SOE"), text.find("$$EOE")
    if soe < 0 or eoe < 0:
        raise ValueError("No $$SOE/$$EOE in response")
    records = []
    for line in text[soe+5:eoe].strip().split("\n"):
        p = [x.strip() for x in line.split(",") if x.strip()]
        if len(p) >= 8:
            try:
                records.append({
                    "jd": float(p[0]), "date": p[1].strip().strip("'\""),
                    "x": float(p[2]), "y": float(p[3]), "z": float(p[4]),
                    "vx": float(p[5]), "vy": float(p[6]), "vz": float(p[7]),
                })
            except ValueError: pass
    return records

# Fetch:

def cmd_fetch(args):
    cfg = read_config()
    step = f"{cfg['output_hours']}h"
    years = cfg["num_years"]

    bodies = PLANETS + (MOONS if args.moons else [])
    stop = (datetime.strptime(args.start, "%Y-%m-%d") +
            timedelta(days=years * 365.25)).strftime("%Y-%m-%d")

    print(f"Config: num_years={years}, output_hours={cfg['output_hours']} -> step={step}")
    print(f"Fetching {len(bodies)} bodies: {args.start} -> {stop}\n")
    data = {}

    for i, (bid, name, gm, parent) in enumerate(bodies):
        print(f"  [{i+1}/{len(bodies)}] {name}...", end=" ", flush=True)
        try:
            recs = parse_vectors(fetch_vectors(bid, args.start, stop, step))
            data[name] = {"id": bid, "gm": gm, "mass": gm/G_KM3,
                          "parent": parent, "states": recs}
            print(f"{len(recs)} epochs")
        except Exception as e:
            print(f"FAILED: {e}")
        time.sleep(0.5)

    if not data:
        print("No data fetched."); return

    # Ensure validation directory exists
    Path("src/validation").mkdir(exist_ok=True)

    with open("src/Body.hpp", "w") as f:
        f.write('#pragma once\n#include "Particle/Particle.hpp"\n#include "Config.hpp"\n\n')
        f.write("struct Body {\n    const char* name;\n    double mass;\n")
        f.write("    double x, y, z;\n    double v_x, v_y, v_z;\n};\n\n")
        f.write(f"// JPL Horizons | {list(data.values())[0]['states'][0]['date']}\n")
        f.write("inline constexpr Body bodies[] = {\n")
        names = list(data.keys())
        for i, n in enumerate(names):
            s = data[n]["states"][0]
            comma = "," if i < len(names)-1 else ""
            f.write(f'    {{ "{n}", {data[n]["mass"]:.6e},\n')
            f.write(f'        {s["x"]:.10e}, {s["y"]:.10e}, {s["z"]:.10e},\n')
            f.write(f'        {s["vx"]:.10e}, {s["vy"]:.10e}, {s["vz"]:.10e} }}{comma}\n')
        f.write("};\n\n")
        f.write("inline void initialize_bodies( Particles &particles, std::size_t const num_bodies ) {\n")
        f.write("    for ( std::size_t i{}; i < num_bodies; ++i ) {\n")
        for attr, field in [("mass","mass"),("pos_x","x"),("pos_y","y"),("pos_z","z"),
                            ("vel_x","v_x"),("vel_y","v_y"),("vel_z","v_z")]:
            mult = " * config::KM_TO_M" if attr != "mass" else ""
            f.write(f'        particles.{attr}()[i] = bodies[i].{field}{mult};\n')
        for a in ["acc_x","acc_y","acc_z"]:
            f.write(f'        particles.{a}()[i] = 0.0;\n')
        f.write("    }\n}\n")
    print(f"\n-> src/Body.hpp ({len(data)} bodies)")

    with open("src/validation/jpl_reference.csv", "w") as f:
        f.write("name,jd,date,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms\n")
        for n, d in data.items():
            for s in d["states"]:
                f.write(f'{n},{s["jd"]:.6f},{s["date"]},'
                        f'{s["x"]:.10e},{s["y"]:.10e},{s["z"]:.10e},'
                        f'{s["vx"]:.10e},{s["vy"]:.10e},{s["vz"]:.10e}\n')
    print(f"-> src/validation/jpl_reference.csv")

    cat = {n: {"id":d["id"],"mass_kg":d["mass"],"gm":d["gm"],"parent":d["parent"],
               "epochs":len(d["states"])} for n,d in data.items()}
    Path("src/validation/body_catalog.json").write_text(json.dumps(cat, indent=2))
    print(f"-> src/validation/body_catalog.json")

# Binary loader:

def load_sim_binary(path):
    data = defaultdict(lambda: {"t": [], "pos": [], "vel": []})
    with open(path, "rb") as f:
        # Read header
        num_bodies = struct.unpack("Q", f.read(8))[0]
        names = []
        for _ in range(num_bodies):
            raw = f.read(32)
            names.append(raw.split(b'\x00')[0].decode())

        # Read frames
        doubles_per_frame = 2 + num_bodies * 6
        frame_size = doubles_per_frame * 8

        while True:
            raw = f.read(frame_size)
            if len(raw) < frame_size:
                break

            frame = np.frombuffer(raw, dtype=np.float64)

            step = struct.unpack("Q", struct.pack("d", frame[0]))[0]
            time_s = frame[1]

            states = frame[2:].reshape(num_bodies, 6)
            for i, name in enumerate(names):
                data[name]["t"].append(time_s)
                data[name]["pos"].append(states[i, 0:3] / 1e3)
                data[name]["vel"].append(states[i, 3:6] / 1e3)

    return {n: {k: np.array(v) for k, v in d.items()} for n, d in data.items()}

# Compare:

def load_ref(path):
    data = defaultdict(lambda: {"jd":[],"pos":[],"vel":[]})
    with open(path) as f:
        for row in csv.DictReader(f):
            n = row["name"]
            data[n]["jd"].append(float(row["jd"]))
            data[n]["pos"].append([float(row["x_km"]),float(row["y_km"]),float(row["z_km"])])
            data[n]["vel"].append([float(row["vx_kms"]),float(row["vy_kms"]),float(row["vz_kms"])])
    return {n: {k: np.array(v) for k,v in d.items()} for n,d in data.items()}

def cmd_compare(args):
    ref = load_ref(args.ref)
    sim = load_sim_binary(args.sim)
    shared = sorted(set(sim) & set(ref))
    if args.bodies:
        shared = [b for b in args.bodies.split(",") if b in shared]
    if not shared:
        print(f"No shared bodies.\n  Sim: {sorted(sim)}\n  Ref: {sorted(ref)}"); return

    print(f"Comparing {len(shared)} bodies\n")

    results = {}
    for name in shared:
        n = min(len(sim[name]["pos"]), len(ref[name]["pos"]))
        if n < 2: continue
        sp = sim[name]["pos"][:n]
        rp = ref[name]["pos"][:n]
        pos_err = np.linalg.norm(sp - rp, axis=1)
        r_ref = np.linalg.norm(rp, axis=1)
        rel = np.where(r_ref > 0, pos_err / r_ref, 0)
        results[name] = rel[1:].max() * 100

    print(f"{'Body':<14} {'Max Relative Error (%)'}")
    print("=" * 40)
    for name in sorted(results, key=results.get, reverse=True):
        print(f"{name:<14} {results[name]:.6f}")

    vals = list(results.values())
    print("=" * 40)
    print(f"{'Worst:':<14} {max(vals):.6f}%")
    print(f"{'Average:':<14} {sum(vals)/len(vals):.6f}%")

# Interface:

def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")

    f = sub.add_parser("fetch")
    f.add_argument("--start", default="1950-01-01")
    f.add_argument("--moons", action="store_true")

    c = sub.add_parser("compare")
    c.add_argument("--sim", default="src/validation/sim_output.bin")
    c.add_argument("--ref", default="src/validation/jpl_reference.csv")
    c.add_argument("--bodies", default=None)

    args = p.parse_args()
    if args.cmd == "fetch": cmd_fetch(args)
    elif args.cmd == "compare": cmd_compare(args)
    else: p.print_help()

if __name__ == "__main__":
    main()