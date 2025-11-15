#!/usr/bin/env python3
"""
Advanced N-Body Initial Conditions Generator
Generates various interesting gravitational scenarios for simulation
"""

import numpy as np
from typing import List, Tuple
import sys

class BodyGenerator:
    def __init__(self, seed=None):
        if seed is not None:
            np.random.seed(seed)
    
    def write_bodies(self, bodies: List[Tuple], filename: str, description: str):
        """Write bodies to CSV file with header"""
        with open(filename, 'w') as f:
            f.write(f"# {description}\n")
            for body in bodies:
                f.write(','.join(map(str, body)) + '\n')
        print(f"Generated {len(bodies)} bodies in '{filename}'")
    
    def binary_system(self, m1=1e30, m2=5e29, separation=1e11, eccentricity=0.3):
        """Generate a binary star system with optional orbiting planets"""
        bodies = []
        
        # Calculate orbital velocities for binary
        G = 6.67430e-11
        total_mass = m1 + m2
        r1 = separation * m2 / total_mass  # Distance of m1 from barycenter
        r2 = separation * m1 / total_mass  # Distance of m2 from barycenter
        
        # Orbital velocity
        v_orb1 = np.sqrt(G * m2 * separation / (r1 * total_mass))
        v_orb2 = np.sqrt(G * m1 * separation / (r2 * total_mass))
        
        # Primary star
        bodies.append((-r1, 0, 0, 0, v_orb1, 0, m1))
        # Secondary star
        bodies.append((r2, 0, 0, 0, -v_orb2, 0, m2))
        
        # Add some planets around the binary
        for i in range(5):
            planet_dist = separation * np.random.uniform(2, 5)
            planet_angle = np.random.uniform(0, 2*np.pi)
            planet_mass = 10**np.random.uniform(24, 27)
            
            px = planet_dist * np.cos(planet_angle)
            py = planet_dist * np.sin(planet_angle)
            pz = np.random.uniform(-separation/10, separation/10)
            
            # Circular orbit velocity around barycenter
            v_planet = np.sqrt(G * total_mass / planet_dist)
            pvx = -v_planet * np.sin(planet_angle)
            pvy = v_planet * np.cos(planet_angle)
            
            bodies.append((px, py, pz, pvx, pvy, 0, planet_mass))
        
        return bodies
    
    def colliding_galaxies(self, n_per_galaxy=50, galaxy_radius=5e12, approach_velocity=50000):
        """Generate two galaxies on a collision course"""
        bodies = []
        
        # Galaxy 1 - centered at negative x, moving right
        for i in range(n_per_galaxy):
            r = galaxy_radius * (np.random.random() ** 0.5)
            theta = np.random.uniform(0, 2*np.pi)
            
            x = -2*galaxy_radius + r * np.cos(theta)
            y = r * np.sin(theta)
            z = np.random.normal(0, galaxy_radius/20)  # Thin disk
            
            # Rotation + approach velocity
            v_rot = 100000 * np.sqrt(r/galaxy_radius)
            vx = approach_velocity - v_rot * np.sin(theta) * 0.5
            vy = v_rot * np.cos(theta)
            vz = np.random.normal(0, 5000)
            
            mass = 10**np.random.uniform(28, 30) if i < 2 else 10**np.random.uniform(24, 27)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # Galaxy 2 - centered at positive x, moving left, tilted
        tilt_angle = np.pi/4  # 45 degree tilt
        for i in range(n_per_galaxy):
            r = galaxy_radius * (np.random.random() ** 0.5)
            theta = np.random.uniform(0, 2*np.pi)
            
            # Create in plane
            x_local = r * np.cos(theta)
            y_local = r * np.sin(theta)
            z_local = np.random.normal(0, galaxy_radius/20)
            
            # Apply tilt around y-axis
            x = 2*galaxy_radius + x_local * np.cos(tilt_angle) - z_local * np.sin(tilt_angle)
            y = y_local
            z = x_local * np.sin(tilt_angle) + z_local * np.cos(tilt_angle)
            
            # Rotation + approach velocity
            v_rot = 100000 * np.sqrt(r/galaxy_radius)
            vx_local = -v_rot * np.sin(theta) * 0.5
            vy_local = v_rot * np.cos(theta)
            vz_local = np.random.normal(0, 5000)
            
            # Apply tilt to velocities
            vx = -approach_velocity + vx_local * np.cos(tilt_angle) - vz_local * np.sin(tilt_angle)
            vy = vy_local
            vz = vx_local * np.sin(tilt_angle) + vz_local * np.cos(tilt_angle)
            
            mass = 10**np.random.uniform(28, 30) if i < 2 else 10**np.random.uniform(24, 27)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies
    
    def planetary_system_with_rings(self, star_mass=2e30):
        """Generate a star with planets, including one with a ring system"""
        bodies = []
        G = 6.67430e-11
        
        # Central star
        bodies.append((0, 0, 0, 0, 0, 0, star_mass))
        
        # Inner rocky planets
        for i in range(3):
            dist = (i + 1) * 5e10
            angle = np.random.uniform(0, 2*np.pi)
            mass = 10**np.random.uniform(24, 25)
            
            x = dist * np.cos(angle)
            y = dist * np.sin(angle)
            z = 0
            
            v = np.sqrt(G * star_mass / dist)
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        # Gas giant with rings (like Saturn)
        giant_dist = 1.5e11
        giant_mass = 5e27
        giant_x = giant_dist
        giant_y = 0
        giant_vx = 0
        giant_vy = np.sqrt(G * star_mass / giant_dist)
        
        bodies.append((giant_x, giant_y, 0, giant_vx, giant_vy, 0, giant_mass))
        
        # Create ring system (many small bodies)
        ring_inner = 7e7
        ring_outer = 1.5e8
        n_ring_particles = 100
        
        for i in range(n_ring_particles):
            r = np.random.uniform(ring_inner, ring_outer)
            theta = np.random.uniform(0, 2*np.pi)
            
            # Position relative to giant
            x = giant_x + r * np.cos(theta)
            y = giant_y + r * np.sin(theta)
            z = np.random.normal(0, 1e6)  # Very thin ring
            
            # Orbital velocity around giant + giant's velocity
            v_ring = np.sqrt(G * giant_mass / r)
            vx = giant_vx - v_ring * np.sin(theta)
            vy = giant_vy + v_ring * np.cos(theta)
            
            mass = 10**np.random.uniform(18, 20)  # Small particles
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        # Outer ice giants
        for i in range(2):
            dist = (3 + i) * 1e11
            angle = np.random.uniform(0, 2*np.pi)
            mass = 10**np.random.uniform(26, 27)
            
            x = dist * np.cos(angle)
            y = dist * np.sin(angle)
            z = np.random.uniform(-1e10, 1e10)
            
            v = np.sqrt(G * star_mass / dist)
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        return bodies
    
    def trojan_asteroids(self, star_mass=2e30, planet_mass=2e27, planet_dist=1.5e11):
        """Generate a system with Trojan asteroids at Lagrange points"""
        bodies = []
        G = 6.67430e-11
        
        # Star
        bodies.append((0, 0, 0, 0, 0, 0, star_mass))
        
        # Planet
        planet_v = np.sqrt(G * star_mass / planet_dist)
        bodies.append((planet_dist, 0, 0, 0, planet_v, 0, planet_mass))
        
        # L4 Trojans (60 degrees ahead)
        l4_angle = np.pi/3
        l4_x = planet_dist * np.cos(l4_angle)
        l4_y = planet_dist * np.sin(l4_angle)
        
        for i in range(30):
            # Scatter around L4 point
            scatter = planet_dist * 0.05
            x = l4_x + np.random.normal(0, scatter)
            y = l4_y + np.random.normal(0, scatter)
            z = np.random.normal(0, scatter/10)
            
            # Velocity similar to planet's orbital velocity
            angle = np.arctan2(y, x)
            v = planet_v * (1 + np.random.normal(0, 0.01))
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            vz = np.random.normal(0, 100)
            
            mass = 10**np.random.uniform(20, 22)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # L5 Trojans (60 degrees behind)
        l5_angle = -np.pi/3
        l5_x = planet_dist * np.cos(l5_angle)
        l5_y = planet_dist * np.sin(l5_angle)
        
        for i in range(30):
            scatter = planet_dist * 0.05
            x = l5_x + np.random.normal(0, scatter)
            y = l5_y + np.random.normal(0, scatter)
            z = np.random.normal(0, scatter/10)
            
            angle = np.arctan2(y, x)
            v = planet_v * (1 + np.random.normal(0, 0.01))
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            vz = np.random.normal(0, 100)
            
            mass = 10**np.random.uniform(20, 22)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies
    
    def hierarchical_triple(self):
        """Generate a hierarchical triple star system"""
        bodies = []
        G = 6.67430e-11
        
        # Inner binary
        m1, m2 = 1.5e30, 1.2e30
        inner_sep = 5e10
        
        # Outer companion
        m3 = 8e29
        outer_sep = 5e11
        
        # Inner binary barycenter is at origin
        r1 = inner_sep * m2 / (m1 + m2)
        r2 = inner_sep * m1 / (m1 + m2)
        
        v1 = np.sqrt(G * m2 * inner_sep / (r1 * (m1 + m2)))
        v2 = np.sqrt(G * m1 * inner_sep / (r2 * (m1 + m2)))
        
        # Inner binary
        bodies.append((-r1, 0, 0, 0, v1, 0, m1))
        bodies.append((r2, 0, 0, 0, -v2, 0, m2))
        
        # Outer companion orbits the inner binary's barycenter
        v3 = np.sqrt(G * (m1 + m2) / outer_sep)
        bodies.append((outer_sep, 0, 0, 0, v3, 0, m3))
        
        # Add some planets in stable orbits
        for i in range(3):
            # Circumbinary planets
            dist = outer_sep * (2 + i * 0.5)
            angle = np.random.uniform(0, 2*np.pi)
            
            x = dist * np.cos(angle)
            y = dist * np.sin(angle)
            z = np.random.uniform(-1e10, 1e10)
            
            v = np.sqrt(G * (m1 + m2 + m3) / dist)
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            
            mass = 10**np.random.uniform(24, 26)
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        return bodies
    
    def globular_cluster(self, n_bodies=1000, radius=1e13, core_radius=2e12):
        """Generate a globular cluster with core collapse dynamics"""
        bodies = []
        G = 6.67430e-11
        total_mass = n_bodies * 1e30  # Approximate
        
        for i in range(n_bodies):
            # King model-like distribution (denser core)
            if np.random.random() < 0.3:  # 30% in core
                r = core_radius * (np.random.random() ** (1/3))
            else:  # 70% in halo
                r = core_radius + (radius - core_radius) * (np.random.random() ** (1/2))
            
            # Random direction
            theta = np.arccos(2 * np.random.random() - 1)
            phi = 2 * np.pi * np.random.random()
            
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            # Velocity dispersion decreases with radius
            sigma = np.sqrt(G * total_mass / (r + core_radius))
            vx = np.random.normal(0, sigma * 0.3)
            vy = np.random.normal(0, sigma * 0.3)
            vz = np.random.normal(0, sigma * 0.3)
            
            # Variable stellar masses (main sequence to giants)
            if i < 5:  # Few massive stars/black holes
                mass = 10**np.random.uniform(30, 31)
            elif i < 20:  # Giants
                mass = 10**np.random.uniform(29, 30)
            else:  # Main sequence
                mass = 10**np.random.uniform(28, 29)
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies
    
    def asteroid_collision(self):
        """Two asteroid fields on collision course"""
        bodies = []
        
        # Field 1 - moving right
        for i in range(40):
            x = np.random.uniform(-2e11, -1e11)
            y = np.random.uniform(-5e10, 5e10)
            z = np.random.uniform(-5e10, 5e10)
            
            vx = np.random.uniform(20000, 30000)
            vy = np.random.uniform(-5000, 5000)
            vz = np.random.uniform(-5000, 5000)
            
            mass = 10**np.random.uniform(20, 23)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        # Field 2 - moving left
        for i in range(40):
            x = np.random.uniform(1e11, 2e11)
            y = np.random.uniform(-5e10, 5e10)
            z = np.random.uniform(-5e10, 5e10)
            
            vx = np.random.uniform(-30000, -20000)
            vy = np.random.uniform(-5000, 5000)
            vz = np.random.uniform(-5000, 5000)
            
            mass = 10**np.random.uniform(20, 23)
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies
    
    def protoplanetary_disk(self, n_particles=150):
        """Young star with protoplanetary disk that will form planets"""
        bodies = []
        G = 6.67430e-11
        star_mass = 2e30
        
        # Young star
        bodies.append((0, 0, 0, 0, 0, 0, star_mass))
        
        # Disk particles with Keplerian rotation
        inner_radius = 3e10
        outer_radius = 5e11
        
        for i in range(n_particles):
            # Power law distribution (more particles closer in)
            r = inner_radius * (outer_radius/inner_radius) ** np.random.random()
            theta = np.random.uniform(0, 2*np.pi)
            
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            # Thin disk with scale height proportional to radius
            z = np.random.normal(0, r * 0.01)
            
            # Keplerian velocity
            v_kep = np.sqrt(G * star_mass / r)
            # Add some radial drift
            v_r = np.random.normal(0, v_kep * 0.01)
            
            vx = -v_kep * np.sin(theta) + v_r * np.cos(theta)
            vy = v_kep * np.cos(theta) + v_r * np.sin(theta)
            vz = np.random.normal(0, v_kep * 0.001)
            
            # Mass distribution (dust to planetesimals)
            if i < 10:  # Few large planetesimals
                mass = 10**np.random.uniform(23, 25)
            else:  # Mostly dust and small particles
                mass = 10**np.random.uniform(18, 22)
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies

def main():
    # Randomly pick a scenario
    gen = BodyGenerator()
    
    scenarios = {
        'binary': (gen.binary_system, "Binary Star System with Planets"),
        'galaxies': (gen.colliding_galaxies, "Colliding Galaxies"),
        'rings': (gen.planetary_system_with_rings, "Planetary System with Ringed Giant"),
        'trojans': (gen.trojan_asteroids, "System with Trojan Asteroids"),
        'triple': (gen.hierarchical_triple, "Hierarchical Triple Star System"),
        'cluster': (gen.globular_cluster, "Globular Cluster"),
        'collision': (gen.asteroid_collision, "Colliding Asteroid Fields"),
        'disk': (gen.protoplanetary_disk, "Protoplanetary Disk")
    }
    
    # Pick a random scenario
    scenario_name = np.random.choice(list(scenarios.keys()))
    func, desc = scenarios[scenario_name]
    
    print(f"Generating scenario: {desc}")
    print("=" * 50)
    
    # Generate the bodies
    bodies = func()
    
    # Write to bodies.csv
    gen.write_bodies(bodies, 'bodies.csv', desc)
    
    # Print some statistics
    masses = [b[6] for b in bodies]
    print(f"\nStatistics:")
    print(f"  Total bodies: {len(bodies)}")
    print(f"  Mass range: {min(masses):.2e} - {max(masses):.2e} kg")
    print(f"  Total mass: {sum(masses):.2e} kg")
    print(f"\nTip: Run again to generate a different random scenario!")

if __name__ == "__main__":
    main()