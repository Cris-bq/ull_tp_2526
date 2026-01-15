# Barnes-Hut N-body Simulation

Barnes-Hut algorithm for efficient gravitational N-body simulations with serial, parallel (OpenMP), and distributed (MPI) implementations.

## Project Structure

```
.
├── Makefile              # Build system
├── README.md             # This file
├── geometry.f90          # 3D vector/point operations
├── particle.f90          # Particle type definition
├── barnes_hut.f90        # Barnes-Hut tree algorithm
├── main.f90              # Serial/OpenMP main program
├── main_mpi.f90          # MPI version
├── comparation_times.py  # Input file generator and comparation times
```

## Algorithm Overview

The Barnes-Hut algorithm reduces N-body force calculation complexity from O(N²) to O(N log N) using an octree:

1. **Tree Construction**: Recursively divide 3D space into octants, creating a tree where each leaf contains ≤1 particle
2. **Mass Calculation**: Calculate total mass and center of mass for each tree node
3. **Force Calculation**: For each particle, traverse the tree:
   - If a node is far enough (l/D < θ), treat it as a single mass
   - Otherwise, recurse into its children
4. **Time Integration**: Update positions and velocities using leapfrog integration

### Key Parameters
- **θ (theta) = 1.0**: Opening angle criterion (smaller = more accurate but slower)
- **Leapfrog integration**: Simple integrator that conserves energy

## Compilation

### Build Commands

```bash
# Compile all versions
make

# Compile specific versions
make serial      # Serial version
make parallel    # OpenMP version
make mpi         # MPI version

# Clean compiled files
make clean

# Show help
make help
```

## Execution

### 1. Serial Version

```bash
./nbody_serial < input.dat
```

### 2. OpenMP Parallel Version

```bash
# Set number of threads
export OMP_NUM_THREADS=4

# Run
./nbody_parallel < input.dat
```

To test with different thread counts:
```bash
export OMP_NUM_THREADS=1 && ./nbody_parallel < input.dat
export OMP_NUM_THREADS=2 && ./nbody_parallel < input.dat
export OMP_NUM_THREADS=4 && ./nbody_parallel < input.dat
export OMP_NUM_THREADS=8 && ./nbody_parallel < input.dat
```

### 3. MPI Distributed Version

```bash
# Run with 4 processes
mpirun --allow-run-as-root -np 4 ./nbody_mpi < input.dat
```

To test with different process counts:
```bash
mpirun --allow-run-as-root -np 1 ./nbody_mpi < input.dat
mpirun --allow-run-as-root -np 2 ./nbody_mpi < input.dat
mpirun --allow-run-as-root -np 4 ./nbody_mpi < input.dat
mpirun --allow-run-as-root -np 8 ./nbody_mpi < input.dat
```

## Creating Input Files

### Input File Format

```
dt          # Time step 
dt_out      # Output interval 
t_end       # End time 
n           # Number of particles 
# For each particle:
mass x y z vx vy vz
```

### Using the Input Generator and benchmark of times

Test all versions with multiple problem sizes:

```bash
python3 comparation_times.py
```

This will:
- Generate input files with 100, 1000, 5000, 10000, 50000 particles
- Run serial, OpenMP (2,4 threads), and MPI (2,4 processes) versions
- Calculate speedups and efficiencies
- Display results table

## Output Files

- `output.dat`: Serial/OpenMP results
- `output_mpi.dat`: MPI results

Format: Each line contains `time x1 y1 z1 x2 y2 z2 ... xn yn zn`

## Parallelization Strategies

### OpenMP (Shared Memory)
- **Parallelized section**: Force calculation loop
- **Strategy**: Each thread calculates forces for a subset of particles
- **Synchronization**: None needed (each particle's acceleration is independent)
- **Best for**: Medium-sized problems on multicore machines

```fortran
!$omp parallel do default(shared) private(i) schedule(dynamic)
do i = 1, n
    call calculate_forces_aux(i, head, particles, accelerations)
end do
!$omp end parallel do
```

### MPI (Distributed Memory)
- **Strategy**: Static domain decomposition
- **Tree**: Each process builds complete tree (needed for force calculations)
- **Communication**: 
  - `MPI_Bcast`: Distribute particle data
  - `MPI_Allreduce`: Sum accelerations from all processes
- **Best for**: Very large problems on clusters

### Performance Characteristics

| Problem Size | Best Choice | Why |
|--------------|-------------|-----|
| N < 500 | Serial or OpenMP | Minimal overhead |
| N = 500-5000 | OpenMP | Good speedup, low overhead |
| N > 5000 | MPI or OpenMP | Depends on hardware |

**Note**: MPI shows overhead for small problems due to communication costs. Performance improves relative to serial as N increases.

## Expected Speedups

Typical results (will vary based on hardware):

| Version | Speedup (N=100) | Speedup (N=1000) |
|---------|-----------------|------------------|
| Serial | 1.0× | 1.0× |
| OpenMP (2T) | 1.7× | 1.8× |
| OpenMP (4T) | 2.5× | 3.2× |
| MPI (2P) | 0.6× | 1.6× |
| MPI (4P) | 0.5× | 2.8× |

## Troubleshooting

### OpenMP not working
- Check compilation: `gfortran --version` should show OpenMP support
- Verify threads: The program prints "OpenMP threads: N" at startup

### MPI errors
- Install OpenMPI: `sudo apt install openmpi-bin libopenmpi-dev`
- Check installation: `mpirun --version`
- Root warning: Use `--allow-run-as-root` flag (safe for development)

### Compilation warnings
- All warnings have been addressed in the latest version
- If you see warnings, run `make clean && make`

## Technical Details

### Barnes-Hut Tree Structure
- **Cell types**: 
  - 0 = empty
  - 1 = leaf (single particle)
  - 2 = internal node (8 children)
- **Tree rebuild**: Complete rebuild every timestep after position updates
- **Memory management**: Proper cleanup with `borrar_tree` to prevent leaks

### Time Integration
Uses leapfrog (kick-drift-kick) method:
```
v(t+dt/2) = v(t) + a(t) × dt/2
r(t+dt) = r(t) + v(t+dt/2) × dt
a(t+dt) = calculate_forces(r(t+dt))
v(t+dt) = v(t+dt/2) + a(t+dt) × dt/2
```

## Authors

Course Exercise 2 - Programming Techniques (Técnicas de Programación)  
Universidad de La Laguna, 2025-2026