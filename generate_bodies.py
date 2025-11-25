#!/usr/bin/env python3
import numpy as np

class BodyGenerator:
    def __init__(self, seed=None):
        if seed is not None:
            np.random.seed(seed)
    
    def generate_globular_cluster(self, n_bodies=1000):
        """Generate a globular cluster with 1000 bodies for benchmarking"""
        bodies = []
        G = 6.67430e-11
        
        # Cluster parameters
        radius = 1e13  # Overall cluster radius
        core_radius = 2e12  # Dense core region
        total_mass = n_bodies * 1e30  # Approximate total mass
        
        print(f"Generating globular cluster with {n_bodies} bodies...")
        print("This configuration is optimized for benchmarking N-body simulations")
        
        for i in range(n_bodies):
            # King model-like distribution (denser core)
            if np.random.random() < 0.3:  # 30% in dense core
                r = core_radius * (np.random.random() ** (1/3))
            else:  # 70% in outer halo
                r = core_radius + (radius - core_radius) * (np.random.random() ** (1/2))
            
            # Random 3D position (spherical distribution)
            theta = np.arccos(2 * np.random.random() - 1)
            phi = 2 * np.pi * np.random.random()
            
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            # Velocity dispersion decreases with radius (virial equilibrium)
            sigma = np.sqrt(G * total_mass / (r + core_radius))
            vx = np.random.normal(0, sigma * 0.3)
            vy = np.random.normal(0, sigma * 0.3)
            vz = np.random.normal(0, sigma * 0.3)
            
            # Mass distribution
            if i < 5:  # 5 supermassive objects (black holes)
                mass = 10**np.random.uniform(30, 31)
            elif i < 20:  # 15 giant stars
                mass = 10**np.random.uniform(29, 30)
            elif i < 100:  # 80 intermediate mass stars
                mass = 10**np.random.uniform(28.5, 29)
            else:  # 900 main sequence stars
                mass = 10**np.random.uniform(28, 28.5)
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies, "Globular Cluster Benchmark - 1000 Bodies"
    
    def generate_planetary_system(self):
        """Generate a star system with multiple planets, rings, and moons in 3D orbits"""
        bodies = []
        G = 6.67430e-11
        
        print("Generating complex planetary system with rings and 3D orbits...")
        print("Enhanced 3D movement with significant Z-axis motion...")
        
        # Central star (Sun-like)
        star_mass = 2e30
        bodies.append((0, 0, 0, 0, 0, 0, star_mass))
        
        # Inner rocky planets with significant inclinations (4 planets)
        rocky_planets = [
            {'dist': 5e10, 'mass': 3e24, 'incl': 0.3, 'z_offset': 0.2},  # High inclination
            {'dist': 1e11, 'mass': 5e24, 'incl': -0.25, 'z_offset': -0.15},  # Retrograde tilt
            {'dist': 1.5e11, 'mass': 6e24, 'incl': 0.4, 'z_offset': 0.3},  # Very tilted
            {'dist': 2e11, 'mass': 4e24, 'incl': -0.35, 'z_offset': -0.25},  # Another retrograde
        ]
        
        for i, planet in enumerate(rocky_planets):
            dist = planet['dist']
            mass = planet['mass']
            inclination = planet['incl']  # Direct inclination value
            z_offset = planet['z_offset']  # Additional Z offset
            angle = i * np.pi / 2 + np.random.uniform(-0.2, 0.2)  # Spread with variation
            
            # 3D position with significant Z component
            x = dist * np.cos(angle) * np.cos(inclination)
            y = dist * np.sin(angle) * np.cos(inclination)
            z = dist * np.sin(inclination) + dist * z_offset * np.sin(angle * 2)
            
            # Orbital velocity with 3D components
            v = np.sqrt(G * star_mass / dist)
            # Add significant Z-component to velocity
            vx = -v * np.sin(angle) * np.cos(inclination)
            vy = v * np.cos(angle) * np.cos(inclination)
            vz = v * np.sin(inclination) * 0.5  # Significant Z velocity
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # Asteroid belt with 3D distribution (15 bodies)
        print("  Adding 3D asteroid belt...")
        for i in range(15):
            dist = np.random.uniform(2.5e11, 3.5e11)
            angle = np.random.uniform(0, 2*np.pi)
            # Much higher inclinations for asteroids
            inclination = np.random.uniform(-0.5, 0.5)  # Up to 45 degrees
            
            x = dist * np.cos(angle) * np.cos(inclination)
            y = dist * np.sin(angle) * np.cos(inclination)
            z = dist * np.sin(inclination)  # Direct Z component
            
            v = np.sqrt(G * star_mass / dist)
            vx = -v * np.sin(angle) * np.cos(inclination)
            vy = v * np.cos(angle) * np.cos(inclination)
            vz = v * np.sin(inclination) * np.random.uniform(0.3, 0.7)  # Variable Z velocity
            
            mass = 10**np.random.uniform(18, 21)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # Gas giant with ring system at an angle (1 giant + 50 ring particles)
        print("  Adding tilted gas giant with 3D ring system...")
        giant_dist = 5e11
        giant_mass = 5.7e27
        giant_incl = 0.4  # 40% tilt
        giant_angle = np.pi / 3
        
        # Giant planet with significant Z position
        giant_x = giant_dist * np.cos(giant_angle)
        giant_y = giant_dist * np.sin(giant_angle) * 0.7
        giant_z = giant_dist * 0.5  # Significant Z offset
        
        v_giant = np.sqrt(G * star_mass / giant_dist)
        giant_vx = -v_giant * np.sin(giant_angle) * 0.8
        giant_vy = v_giant * np.cos(giant_angle) * 0.8
        giant_vz = v_giant * 0.3  # Z velocity component
        
        bodies.append((giant_x, giant_y, giant_z, giant_vx, giant_vy, giant_vz, giant_mass))
        
        # Ring system tilted in 3D (50 particles)
        for i in range(50):
            r = np.random.uniform(7e7, 1.5e8)
            theta = np.random.uniform(0, 2*np.pi)
            
            # Ring tilted at multiple angles for 3D effect
            ring_tilt_x = 0.4  # Tilt in X
            ring_tilt_y = 0.3  # Tilt in Y
            ring_tilt_z = 0.5  # Tilt in Z
            
            # 3D rotation of ring particles
            local_x = r * np.cos(theta)
            local_y = r * np.sin(theta) * np.cos(ring_tilt_x)
            local_z = r * np.sin(theta) * np.sin(ring_tilt_y) + r * 0.2 * np.cos(theta * 2)
            
            # Transform to global coordinates
            x = giant_x + local_x * np.cos(ring_tilt_z) - local_y * np.sin(ring_tilt_z)
            y = giant_y + local_x * np.sin(ring_tilt_z) + local_y * np.cos(ring_tilt_z)
            z = giant_z + local_z + np.random.normal(0, 2e6)  # Some thickness
            
            # Orbital velocity around giant + giant's velocity
            v_ring = np.sqrt(G * giant_mass / r)
            vx = giant_vx - v_ring * np.sin(theta) * np.cos(ring_tilt_x)
            vy = giant_vy + v_ring * np.cos(theta) * np.cos(ring_tilt_y)
            vz = giant_vz + v_ring * 0.2 * np.sin(theta * 2)  # Complex Z motion
            
            mass = 10**np.random.uniform(16, 19)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # Moons with inclined orbits (3 moons)
        moon_params = [
            {'dist': 2e8, 'incl': 0.3},
            {'dist': 4e8, 'incl': -0.4},
            {'dist': 6e8, 'incl': 0.5},
        ]
        
        for i, params in enumerate(moon_params):
            moon_dist = params['dist']
            incl = params['incl']
            moon_angle = i * 2 * np.pi / 3
            
            # 3D moon positions
            moon_x = giant_x + moon_dist * np.cos(moon_angle) * np.cos(incl)
            moon_y = giant_y + moon_dist * np.sin(moon_angle) * np.cos(incl)
            moon_z = giant_z + moon_dist * np.sin(incl)
            
            v_moon = np.sqrt(G * giant_mass / moon_dist)
            moon_vx = giant_vx - v_moon * np.sin(moon_angle) * np.cos(incl)
            moon_vy = giant_vy + v_moon * np.cos(moon_angle) * np.cos(incl)
            moon_vz = giant_vz + v_moon * np.sin(incl) * 0.5
            
            moon_mass = 10**np.random.uniform(22, 23)
            bodies.append((moon_x, moon_y, moon_z, moon_vx, moon_vy, moon_vz, moon_mass))
        
        # Second gas giant with polar orbit (1 body)
        jupiter_dist = 7.8e11
        jupiter_mass = 1.9e27
        # Polar orbit - moving mostly in XZ plane
        jupiter_x = jupiter_dist * 0.3
        jupiter_y = jupiter_dist * 0.2
        jupiter_z = jupiter_dist * 0.8  # High Z component
        
        v_jupiter = np.sqrt(G * star_mass / jupiter_dist)
        jupiter_vx = v_jupiter * 0.2
        jupiter_vy = v_jupiter * 0.3
        jupiter_vz = -v_jupiter * 0.7  # Strong Z velocity (polar orbit)
        
        bodies.append((jupiter_x, jupiter_y, jupiter_z, jupiter_vx, jupiter_vy, jupiter_vz, jupiter_mass))
        
        # Ice giants with highly inclined orbits (2 bodies)
        ice_giants = [
            {'dist': 1.5e12, 'mass': 8.7e25, 'polar': True},  # Polar orbit
            {'dist': 2e12, 'mass': 1.0e26, 'polar': False},  # Diagonal orbit
        ]
        
        for i, giant in enumerate(ice_giants):
            dist = giant['dist']
            mass = giant['mass']
            
            if giant['polar']:
                # Polar orbit (XZ plane)
                angle = np.pi * i
                x = dist * np.cos(angle) * 0.2
                y = dist * 0.1
                z = dist * np.sin(angle) * 0.9
                
                v = np.sqrt(G * star_mass / dist)
                vx = -v * np.sin(angle) * 0.2
                vy = v * 0.1
                vz = v * np.cos(angle) * 0.9
            else:
                # Diagonal orbit (equal components)
                x = dist / np.sqrt(3)
                y = dist / np.sqrt(3)
                z = dist / np.sqrt(3)
                
                v = np.sqrt(G * star_mass / dist)
                vx = -v / np.sqrt(3)
                vy = v / np.sqrt(3)
                vz = v / np.sqrt(3)
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # Comet-like objects with extreme 3D orbits (22 bodies to total 100)
        print("  Adding comet-like objects with extreme 3D orbits...")
        for i in range(22):
            dist = np.random.uniform(8e11, 3e12)
            
            # Random 3D position on sphere
            theta = np.arccos(2 * np.random.random() - 1)  # Full sphere coverage
            phi = 2 * np.pi * np.random.random()
            
            x = dist * np.sin(theta) * np.cos(phi)
            y = dist * np.sin(theta) * np.sin(phi)
            z = dist * np.cos(theta)  # Full Z range
            
            # Velocity perpendicular to position for orbit
            v = np.sqrt(G * star_mass / dist) * np.random.uniform(0.5, 1.2)
            
            # Random orbit direction for variety
            orbit_theta = np.random.uniform(0, np.pi)
            orbit_phi = np.random.uniform(0, 2*np.pi)
            
            vx = v * np.sin(orbit_theta) * np.cos(orbit_phi)
            vy = v * np.sin(orbit_theta) * np.sin(orbit_phi)
            vz = v * np.cos(orbit_theta)
            
            mass = 10**np.random.uniform(19, 22)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies, "Complex 3D Planetary System (100 bodies)"
    
    def write_bodies(self, bodies, description, filename='bodies.csv'):
        """Write bodies to CSV file"""
        with open(filename, 'w') as f:
            f.write(f"# {description}\n")
            f.write("# Format: x, y, z, vx, vy, vz, mass\n")
            for body in bodies:
                f.write(','.join(map(str, body)) + '\n')
        
        # Statistics
        masses = [b[6] for b in bodies]
        print(f"\nStatistics:")
        print(f"  Total bodies: {len(bodies)}")
        print(f"  Mass range: {min(masses):.2e} - {max(masses):.2e} kg")
        
        if "Planetary" in description:
            print(f"  Star mass: {bodies[0][6]:.2e} kg")
            print("\nObject breakdown (100 total bodies):")
            print(f"  1 star")
            print(f"  4 rocky planets")
            print(f"  15 asteroids")
            print(f"  1 gas giant with 50 ring particles")
            print(f"  3 moons")
            print(f"  1 Jupiter-like giant")
            print(f"  2 ice giants")
            print(f"  22 comet-like objects")
            print("\nNote: Enhanced 3D movement with polar and inclined orbits")
        else:
            print(f"  Total mass: {sum(masses):.2e} kg")
            print(f"  Average mass: {np.mean(masses):.2e} kg")
            
            # Position statistics
            positions = np.array([[b[0], b[1], b[2]] for b in bodies])
            distances = np.sqrt(np.sum(positions**2, axis=1))
            print(f"  Distance range: {min(distances):.2e} - {max(distances):.2e} m")
            print(f"  Mean distance: {np.mean(distances):.2e} m")
        
        print(f"\nFile '{filename}' created successfully!")

def main():
    generator = BodyGenerator()
    
    # Randomly choose between the two configurations
    choice = np.random.choice(['benchmark', 'planetary'])
    
    print("=" * 60)
    print("N-BODY SYSTEM GENERATOR")
    print("=" * 60)
    print(f"\nRandomly selected: {choice.upper()}")
    print("=" * 60)
    
    # Generate based on random choice
    if choice == 'benchmark':
        bodies, description = generator.generate_globular_cluster(n_bodies=1000)
    else:  # planetary
        bodies, description = generator.generate_planetary_system()
    
    # Write to file
    generator.write_bodies(bodies, description)
    
    print("\nRun again to generate a different random configuration!")

if __name__ == "__main__":
    main()