# Simulation of Spinodal Decomposition in Binary Alloys Using Vacancy-Mediated Dynamics

This project simulates spinodal decomposition in binary alloys using vacancy-mediated dynamics, based on the Ising model. Spinodal decomposition occurs when a rapidly quenched binary alloy separates into distinct domains of its two constituent metals. This simulation provides insights into the growth of these domains over time and compares the vacancy-mediated dynamics to spin-exchange dynamics.

## Key Steps in the Project

1. **Initial Configuration**:  
   The simulation begins by placing A and B atoms randomly on a 2D square lattice (128x128) with a single vacancy. This configuration represents a high-energy, infinite-temperature state where the atoms are completely mixed. The vacancy will later drive the dynamics of the system by moving atoms and forming domains.

2. **Quenching**:  
   The system is instantaneously quenched to a low temperature, below the critical temperature \( T_c \), using the Metropolis algorithm. This sudden temperature drop causes the system to leave the mixed state and form separate domains of A and B atoms.

3. **Domain Growth Measurement**:  
   To analyze how domains grow over time, two methods are used:
   - **Pair correlation function** \( C(r) \): Measures the spatial correlation between A and B atoms, with the first zero-crossing of \( C(r) \) corresponding to the domain size.
   - **Energy-based measure** \( R \): Calculated using the system’s energy and perimeter, providing another estimate of domain size.

4. **Scaling Analysis**:  
   The growth of the domain size is expected to follow the law \( R \propto t^{1/3} \). A log-log plot of domain size \( R \) versus time \( t \) will be used to verify this scaling, and the simulation will be repeated at different temperatures (e.g., \( 0.2T_c \), \( 0.7T_c \)) to observe variations in behavior.

5. **Vacancy-Mediated vs Spin-Exchange Dynamics**:  
   The primary focus of this project is vacancy-mediated dynamics, where a single vacancy exchanges places with neighboring atoms to drive domain formation. This approach is more realistic for modeling real-world alloys compared to spin-exchange dynamics, which involve direct atom swapping. 

## Running the Simulation

### File Descriptions

- **`binary_alloys.cpp`**:  
   This program simulates 2D vacancy-mediated dynamics in a binary alloy. It outputs configuration changes over time into text files, as well as calculating and storing values for the pair correlation function \( C(r) \) and domain size \( R \). You can modify parameters such as lattice size and the ratio of A to B atoms directly in the code.

- **`binary_alloys_3d.cpp`**:  
   A 3D version of the above simulation, simulating the same vacancy-mediated dynamics on a 3D lattice. The output format is similar, with results stored in text files.

- **`plot_correlation_func.py`**:  
   This script generates plots of the pair correlation function \( C(r) \), showing how the correlation between atoms changes over time. This helps to track domain formation in the system.

- **`plot_lattice.py`**:  
   This script visualizes the configuration of the lattice at different time steps, providing a clear view of how the domains of A and B atoms evolve over time. *Note: This script only works for 2D simulations.*

- **`plot_R_from_Energy.py`**:  
   This script plots the domain size \( R \) versus time \( t \), based on the energy calculations from the simulation. 
- **`plot_R_from_pair_cor.py`**:  
   This script also plots \( R \) versus \( t \), but uses the pair correlation function \( C(r) \) to determine domain size. Both energy-based and correlation-based measures are compared to validate results.

### Usage Instructions

1. **Compiling and Running**:  
   Compile the C++ code using a standard C++ compiler (e.g., `g++`).
   ```bash
   g++ binary_alloys.cpp -o binary_alloys
   ./binary_alloys
