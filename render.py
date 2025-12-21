from vpython import *
import pandas as pd
import numpy as np
import random

data = pd.read_csv('trajectories.csv')
cols = ['x', 'y', 'z', 'body_id', 'step']
for col in cols:
    if col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
data.dropna(inplace=True)
data['body_id'] = data['body_id'].astype(int)
data['step'] = data['step'].astype(int)
body_ids = data['body_id'].unique()

bodies_data = pd.read_csv('bodies.csv', comment='#', header=None,
                          names=['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'])
body_masses = bodies_data['mass'].values

print(f"Loaded {len(data)} trajectory points for {len(body_ids)} bodies.")
print(f"Loaded masses for {len(body_masses)} bodies from bodies.csv")

all_x = data['x'].values
all_y = data['y'].values
all_z = data['z'].values

max_extent = max(
    np.abs(all_x).max(),
    np.abs(all_y).max(),
    np.abs(all_z).max()
)

scene_range = max_extent * 1.3

first_step = data[data['step'] == data['step'].min()]
positions = []
for b_id in body_ids:
    row = first_step[first_step['body_id'] == b_id]
    if not row.empty:
        positions.append(np.array([row['x'].values[0], row['y'].values[0], row['z'].values[0]]))

min_separation = float('inf')
for i in range(len(positions)):
    for j in range(i + 1, len(positions)):
        dist = np.linalg.norm(positions[i] - positions[j])
        if dist > 0:
            min_separation = min(min_separation, dist)

min_radius = min_separation * 0.1
max_radius = min_separation * 0.5

print(f"Max extent: {max_extent:.2e} m")
print(f"Min separation: {min_separation:.2e} m")
print(f"Radius range: {min_radius:.2e} - {max_radius:.2e} m")

masses = body_masses[body_ids]
min_mass = masses.min()
max_mass = masses.max()

def calculate_radius(mass):
    if max_mass == min_mass:
        return (min_radius + max_radius) / 2
    log_mass = np.log10(mass)
    log_min = np.log10(min_mass)
    log_max = np.log10(max_mass)
    normalized = (log_mass - log_min) / (log_max - log_min)
    return min_radius + normalized * (max_radius - min_radius)

body_radii = [calculate_radius(masses[i]) for i in range(len(body_ids))]

scene.title = "<h2 style='color: #aaa; font-family: Helvetica;'>N-Body Gravitational Simulation</h2>"
scene.background = color.black
scene.width = 1600
scene.height = 900
scene.resizable = True
scene.autoscale = False
scene.range = scene_range
scene.userzoom = True
scene.userspin = True
scene.userpan = True
scene.center = vector(0, 0, 0)

scene.caption = """
<style>
    body { background: #0a0a0f; }
    .controls { 
        font-family: 'Segoe UI', Arial, sans-serif; 
        color: #888; 
        padding: 15px;
        line-height: 1.8;
    }
    .key { 
        background: #222; 
        padding: 3px 8px; 
        border-radius: 4px; 
        color: #fff;
        font-family: monospace;
    }
</style>
<div class="controls">
    <span class="key">SPACE</span> Pause/Resume &nbsp;&nbsp;
    <span class="key">↑↓</span> Speed &nbsp;&nbsp;
    <span class="key">R</span> Reset View &nbsp;&nbsp;
    <span class="key">S</span> Restart Simulation
</div>
"""

scene.lights = []
scene.ambient = color.gray(0.03)
distant_light(direction=vector(0.5, 0.3, 1), color=color.gray(0.4))
distant_light(direction=vector(-0.3, -0.5, -0.5), color=color.gray(0.15))

def create_starfield(num_stars=800, spread=None):
    if spread is None:
        spread = scene_range * 100
    
    stars = []
    for _ in range(num_stars):
        distance = random.uniform(spread * 0.3, spread)
        theta = random.uniform(0, 2 * np.pi)
        phi = random.uniform(-np.pi/2, np.pi/2)
        
        pos = vector(
            distance * np.cos(phi) * np.cos(theta),
            distance * np.cos(phi) * np.sin(theta),
            distance * np.sin(phi)
        )
        
        brightness = random.uniform(0.15, 0.9)
        size = scene_range * random.uniform(0.001, 0.004)
        
        temp_variation = random.uniform(-0.1, 0.1)
        star_color = vector(
            min(1, 0.9 + temp_variation),
            min(1, 0.9 + temp_variation * 0.5),
            min(1, 1.0 - temp_variation * 0.3)
        )
        
        star = sphere(
            pos=pos,
            radius=size,
            color=star_color,
            emissive=True,
            opacity=brightness
        )
        stars.append(star)
    return stars

print("Creating starfield...")
background_stars = create_starfield()

astro_palette = [
    vector(1.0, 0.95, 0.6),    # Sun
    vector(0.75, 0.7, 0.65),   # Mercury
    vector(0.95, 0.85, 0.6),   # Venus 
    vector(0.25, 0.5, 0.95),   # Earth
    vector(0.9, 0.45, 0.25),   # Mars
    vector(0.85, 0.75, 0.55),  # Jupiter
    vector(0.9, 0.82, 0.55),   # Saturn
    vector(0.6, 0.85, 0.92),   # Uranus
    vector(0.35, 0.45, 0.85),  # Neptune
    vector(0.8, 0.75, 0.7),    # Pluto
]

def get_trail_color(base_color):
    return base_color * 0.35

spheres = []
glows = []
star_idx = np.argmax(body_masses[body_ids])

def create_bodies():
    body_spheres = []
    body_glows = []
    
    for idx, b_id in enumerate(body_ids):
        col = astro_palette[idx % len(astro_palette)]
        is_star = (idx == star_idx)
        radius = body_radii[idx]
        
        if is_star:
            core = sphere(
                radius=radius,
                color=col,
                emissive=True,
                shininess=0
            )
            body_spheres.append(core)
            
            glow1 = sphere(
                radius=radius * 1.3,
                color=col * 0.9,
                opacity=0.25,
                emissive=True
            )
            glow2 = sphere(
                radius=radius * 1.8,
                color=col * 0.7,
                opacity=0.1,
                emissive=True
            )
            glow3 = sphere(
                radius=radius * 2.5,
                color=vector(1, 0.9, 0.7) * 0.5,
                opacity=0.04,
                emissive=True
            )
            body_glows.append((idx, [glow1, glow2, glow3]))
        else:
            obj = sphere(
                radius=radius,
                color=col,
                emissive=False,
                shininess=0.4,
                make_trail=True,
                trail_type="curve",
                trail_color=get_trail_color(col),
                trail_radius=radius * 0.08,
                retain=200,
                interval=2
            )
            body_spheres.append(obj)
    
    return body_spheres, body_glows

def destroy_bodies():
    global spheres, glows
    
    for obj in spheres:
        obj.visible = False
        if hasattr(obj, 'clear_trail'):
            obj.clear_trail()
        del obj
    
    for _, glow_list in glows:
        for g in glow_list:
            g.visible = False
            del g
    
    spheres = []
    glows = []

spheres, glows = create_bodies()

ecliptic = ring(
    pos=vector(0, 0, 0),
    axis=vector(0, 0, 1),
    radius=scene_range * 0.75,
    thickness=scene_range * 0.002,
    color=vector(0.3, 0.4, 0.5),
    opacity=0.08
)

steps = sorted(data['step'].unique())
max_step = max(steps)

info_label = label(
    pos=vector(0, 0, 0),
    pixel_pos=True,
    xoffset=-scene.width/2 + 20,
    yoffset=scene.height/2 - 25,
    text='',
    color=color.white,
    opacity=0.6,
    box=False,
    font='monospace',
    height=14,
    align='left'
)

running = True
animation_speed = 30
current_step = 0
restart_requested = False

def keydown(evt):
    global running, animation_speed, restart_requested
    
    if evt.key == ' ':
        running = not running
    elif evt.key == 'up':
        animation_speed = min(500, animation_speed + 20)
    elif evt.key == 'down':
        animation_speed = max(10, animation_speed - 20)
    elif evt.key == 'r':
        scene.center = vector(0, 0, 0)
        scene.range = scene_range
        scene.forward = vector(0, 0, -1)
        scene.up = vector(0, 1, 0)
    elif evt.key == 's':
        restart_requested = True

scene.bind('keydown', keydown)

def reset_simulation():
    global current_step, spheres, glows, running, restart_requested
    
    destroy_bodies()
    spheres, glows = create_bodies()
    current_step = 0
    running = True
    restart_requested = False
    
    info_label.color = color.white
    print("Simulation restarted.")

def update_positions(step_index):
    s = steps[step_index]
    step_data = data[data['step'] == s]
    
    for i, b_id in enumerate(body_ids):
        row = step_data[step_data['body_id'] == b_id]
        if not row.empty:
            new_pos = vector(
                row['x'].values[0],
                row['y'].values[0],
                row['z'].values[0]
            )
            spheres[i].pos = new_pos
            
            for glow_idx, glow_list in glows:
                if glow_idx == i:
                    for g in glow_list:
                        g.pos = new_pos

def update_info_display():
    s = steps[current_step] if current_step < len(steps) else max_step
    progress = (current_step / len(steps)) * 100
    status = "▶ RUNNING" if running else "⏸ PAUSED"
    info_label.text = f"{status}  |  Step: {s}/{max_step}  |  Progress: {progress:.1f}%  |  Speed: {animation_speed}"

print("Starting simulation...")
print("Press 'S' to restart simulation at any time.")

while True:
    rate(animation_speed)
    
    if restart_requested:
        reset_simulation()
        continue
    
    if current_step < len(steps):
        if running:
            update_positions(current_step)
            update_info_display()
            current_step += 1
        else:
            update_info_display()
    else:
        info_label.text = "✓ COMPLETE  |  Press 'S' to restart  |  Speed: " + str(animation_speed)
        info_label.color = vector(0.5, 1, 0.5)