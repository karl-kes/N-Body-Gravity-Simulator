from vpython import *
import pandas as pd
import numpy as np

# 1. Load trajectory data
data = pd.read_csv('trajectories.csv')
cols = ['x', 'y', 'z', 'body_id', 'step']
for col in cols:
    if col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
data.dropna(inplace=True)
data['body_id'] = data['body_id'].astype(int)
data['step'] = data['step'].astype(int)
body_ids = data['body_id'].unique()

# 2. Load masses
bodies_data = pd.read_csv('bodies.csv', comment='#', header=None,
                          names=['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'])
body_masses = bodies_data['mass'].values
print(f"Loaded {len(data)} trajectory points for {len(body_ids)} bodies.")
print(f"Loaded masses for {len(body_masses)} bodies from bodies.csv")

# 3. Calculate sphere radii based on mass
masses = body_masses[body_ids]
min_mass = masses.min()
max_mass = masses.max()
min_radius = 2e8
max_radius = 2e10

def calculate_radius(mass):
    if max_mass == min_mass:
        return (min_radius + max_radius) / 2
    
    log_mass = np.log10(mass)
    log_min = np.log10(min_mass)
    log_max = np.log10(max_mass)
    
    normalized = (log_mass - log_min) / (log_max - log_min)
    radius = min_radius + normalized * (max_radius - min_radius)
    return radius

body_radii = [calculate_radius(masses[i]) for i in range(len(body_ids))]

# 4. Setup the Scene
scene.background = color.black
scene.width = 1400
scene.height = 700
scene.resizable = True
scene.autoscale = False
scene.range = 5e11
scene.userzoom = True
scene.userspin = True
scene.userpan = True

# 5. Create Coordinate Axes
axis_length = -scene.range * 1e3
axis_thickness = scene.range * 0.002
axis_opacity = 0.6

x_axis = arrow(pos=vector(0, 0, 0), 
               axis=vector(axis_length, 0, 0),
               shaftwidth=axis_thickness,
               color=color.red,
               opacity=axis_opacity)
y_axis = arrow(pos=vector(0, 0, 0),
               axis=vector(0, axis_length, 0),
               shaftwidth=axis_thickness,
               color=color.green,
               opacity=axis_opacity)
z_axis = arrow(pos=vector(0, 0, 0),
               axis=vector(0, 0, axis_length),
               shaftwidth=axis_thickness,
               color=color.blue,
               opacity=axis_opacity)

# 6. Create Objects
spheres = []
colors_palette = [color.red, color.cyan, color.yellow, color.green, 
                  color.magenta, color.orange, color.white, color.purple]
for idx, b_id in enumerate(body_ids):
    col = colors_palette[idx % len(colors_palette)]
    obj = sphere(
        radius=body_radii[idx],
        color=col, 
        make_trail=True, 
        trail_type="curve", 
        retain=10
    )
    spheres.append(obj)

# 7. Animation Loop
running = True
animation_speed = 90
steps = data['step'].unique()
current_step = 0

def keydown(evt):
    global running
    if evt.key == ' ':
        running = not running

scene.bind('keydown', keydown)

# Main animation loop - runs only once through all steps
while current_step < len(steps):
    rate(animation_speed)
    
    if running:
        s = steps[current_step]
        step_data = data[data['step'] == s]
        
        for i, b_id in enumerate(body_ids):
            row = step_data[step_data['body_id'] == b_id]
            if not row.empty:
                new_pos = vector(row['x'].values[0], row['y'].values[0], row['z'].values[0])
                spheres[i].pos = new_pos
        current_step += 1

print("Simulation finished.")
while True:
    rate(30)