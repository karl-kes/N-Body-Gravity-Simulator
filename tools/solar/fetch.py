"""
Solar system initial-condition generator (JPL Horizons fetcher).

Queries the NASA JPL Horizons API for all 35 bodies (Sun, 8 planets, Pluto,
25 natural satellites), writes them to tests/sim_ic.bin in the unified binary
IC format that main.cpp reads at runtime. Also writes a reference ephemeris
CSV used by tools/solar/compare.py.

Usage:
    python tools/solar/fetch.py
    python tools/solar/fetch.py --moons
    python tools/solar/fetch.py --start 1900-01-01

Binary IC format (matches main.cpp and tools/galaxy/generate.py):
    Header:
        uint64_t N
        double   dt          (s)
        double   total_time  (s)
        double   output_dt   (s)
        double   eps         (m)
    Names (32 * N bytes, null-padded ASCII)
    Particle data: 7 doubles per particle (mass, x, y, z, vx, vy, vz)  SI units
"""

import argparse, json, ssl, struct, time, urllib.parse, urllib.request
from datetime import datetime, timedelta
from pathlib import Path

KM_TO_M    = 1e3
SEC_PER_YR = 365.25 * 86400.0
SEC_PER_HR = 3600.0
G_KM3      = 6.6743e-20    # km^3 kg^-1 s^-2

# Solar scenario parameters:
SOLAR_DT          = 900.0                       # 15 min (s)
SOLAR_NUM_YEARS   = 249
SOLAR_OUTPUT_HRS  = 487
SOLAR_EPS         = 1e-9
SOLAR_TOTAL_TIME  = SOLAR_NUM_YEARS * SEC_PER_YR
SOLAR_OUTPUT_DT   = SOLAR_OUTPUT_HRS * SEC_PER_HR

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

API = "https://ssd.jpl.nasa.gov/api/horizons.api"

# Unverified SSL for environments where the JPL chain isn't trusted (university/corp).
_ssl_ctx = ssl.create_default_context()
_ssl_ctx.check_hostname = False
_ssl_ctx.verify_mode = ssl.CERT_NONE


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
            with urllib.request.urlopen(url, timeout=60, context=_ssl_ctx) as r:
                return r.read().decode()
        except Exception as e:
            if attempt < 2: time.sleep(2 * (attempt + 1))
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
                    "x": float(p[2]),  "y": float(p[3]),  "z": float(p[4]),
                    "vx": float(p[5]), "vy": float(p[6]), "vz": float(p[7]),
                })
            except ValueError: pass
    return records


def write_ic_file(path, data, names, dt, total_time, output_dt, eps):
    """Write the unified binary IC format."""
    N = len(names)
    Path(path).parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f:
        # Header
        f.write(struct.pack("Q", N))
        f.write(struct.pack("d", dt))
        f.write(struct.pack("d", total_time))
        f.write(struct.pack("d", output_dt))
        f.write(struct.pack("d", eps))

        # Names (32 bytes each, null-padded)
        for name in names:
            buf = name.encode("ascii", errors="replace")[:31]
            f.write(buf + b"\x00" * (32 - len(buf)))

        # Particle data: mass, x, y, z, vx, vy, vz in SI
        for name in names:
            d = data[name]
            s = d["states"][0]
            f.write(struct.pack("d", d["mass"]))
            f.write(struct.pack("d", s["x"]  * KM_TO_M))
            f.write(struct.pack("d", s["y"]  * KM_TO_M))
            f.write(struct.pack("d", s["z"]  * KM_TO_M))
            f.write(struct.pack("d", s["vx"] * KM_TO_M))
            f.write(struct.pack("d", s["vy"] * KM_TO_M))
            f.write(struct.pack("d", s["vz"] * KM_TO_M))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start",  default="1950-01-01")
    p.add_argument("--moons",  action="store_true")
    p.add_argument("--output", default="tests/sim_ic.bin")
    args = p.parse_args()

    step = f"{SOLAR_OUTPUT_HRS}h"
    bodies = PLANETS + (MOONS if args.moons else [])
    stop = (datetime.strptime(args.start, "%Y-%m-%d") +
            timedelta(days=SOLAR_NUM_YEARS * 365.25)).strftime("%Y-%m-%d")

    print(f"Scenario: solar system ({len(bodies)} bodies, {SOLAR_NUM_YEARS} years, dt={SOLAR_DT}s)")
    print(f"Fetching {args.start} -> {stop}\n")
    data = {}

    for i, (bid, name, gm, parent) in enumerate(bodies):
        print(f"  [{i+1}/{len(bodies)}] {name}...", end=" ", flush=True)
        try:
            recs = parse_vectors(fetch_vectors(bid, args.start, stop, step))
            data[name] = {"id": bid, "gm": gm, "mass": gm / G_KM3,
                          "parent": parent, "states": recs}
            print(f"{len(recs)} epochs")
        except Exception as e:
            print(f"FAILED: {e}")
        time.sleep(0.5)

    if not data:
        print("No data fetched.")
        return

    Path("tests").mkdir(exist_ok=True)
    names = list(data.keys())

    write_ic_file(args.output, data, names,
                  SOLAR_DT, SOLAR_TOTAL_TIME, SOLAR_OUTPUT_DT, SOLAR_EPS)
    print(f"\n-> {args.output} ({len(names)} bodies)")

    # Reference CSV consumed by tools/solar/compare.py:
    with open("tests/jpl_reference.csv", "w") as f:
        f.write("name,jd,date,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms\n")
        for n, d in data.items():
            for s in d["states"]:
                f.write(f'{n},{s["jd"]:.6f},{s["date"]},'
                        f'{s["x"]:.10e},{s["y"]:.10e},{s["z"]:.10e},'
                        f'{s["vx"]:.10e},{s["vy"]:.10e},{s["vz"]:.10e}\n')
    print(f"-> tests/jpl_reference.csv")

    cat = {n: {"id": d["id"], "mass_kg": d["mass"], "gm": d["gm"], "parent": d["parent"],
               "epochs": len(d["states"])} for n, d in data.items()}
    Path("tests/body_catalog.json").write_text(json.dumps(cat, indent=2))
    print(f"-> tests/body_catalog.json")


if __name__ == "__main__":
    main()
