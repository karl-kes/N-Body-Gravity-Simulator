from vpython import *
import pandas as pd

# 1. Load Data
data = pd.read_csv('trajectories.csv')
cols = ['x', 'y', 'z', 'body_id', 'step']
for col in cols:
    if col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
data.dropna(inplace=True)
data['body_id'] = data['body_id'].astype(int)
data['step'] = data['step'].astype(int)
body_ids = data['body_id'].unique()
print(f"Loaded {len(data)} valid rows for {len(body_ids)} bodies.")

# 2. Setup the Scene with Resizable Window
scene.background = color.black
scene.width = 1400  # Initial width
scene.height = 800  # Initial height
scene.resizable = True  # Enable window resizing
scene.autoscale = False  # Disable autoscaling to maintain view
scene.range = 5e11
scene.userzoom = True  # Allow user to zoom
scene.userspin = True  # Allow user to rotate view
scene.userpan = True  # Allow user to pan view

# 3. Create Coordinate Axes
axis_length = -scene.range * 10
axis_thickness = scene.range * 0.002
axis_opacity = 0.6

# Function to calculate label size based on scene
def get_label_size():
    return int(scene.width / 50)  # Dynamic label size based on window width

# X-axis (Red)
x_axis = arrow(pos=vector(0, 0, 0), 
               axis=vector(axis_length, 0, 0),
               shaftwidth=axis_thickness,
               color=color.red,
               opacity=axis_opacity)
x_label = label(pos=vector(axis_length, 0, 0),
                text='X',
                color=color.red,
                opacity=0.8,
                height=get_label_size(),
                box=False)

# Y-axis (Green)
y_axis = arrow(pos=vector(0, 0, 0),
               axis=vector(0, axis_length, 0),
               shaftwidth=axis_thickness,
               color=color.green,
               opacity=axis_opacity)
y_label = label(pos=vector(0, axis_length, 0),
                text='Y',
                color=color.green,
                opacity=0.8,
                height=get_label_size(),
                box=False)

# Z-axis (Blue)
z_axis = arrow(pos=vector(0, 0, 0),
               axis=vector(0, 0, axis_length),
               shaftwidth=axis_thickness,
               color=color.blue,
               opacity=axis_opacity)
z_label = label(pos=vector(0, 0, axis_length),
                text='Z',
                color=color.blue,
                opacity=0.8,
                height=get_label_size(),
                box=False)

# 4. Create Objects
spheres = []
colors = [color.red, color.cyan, color.yellow, color.green, color.magenta, color.orange, color.white]
for idx, b_id in enumerate(body_ids):
    col = colors[idx % len(colors)]
    obj = sphere(
        radius=6e8,
        color=col, 
        make_trail=True, 
        trail_type="curve", 
        retain=100
    )
    spheres.append(obj)

# 5. Animation Controls
running = True
animation_speed = 120  # Default frame rate

# 6. Animation Loop
steps = data['step'].unique()
current_step = 0

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
        
        # Loop back to beginning when finished
        if current_step >= len(steps):
            current_step = 0
            print("Simulation restarting...")

print("Simulation finished.")