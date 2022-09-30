# 1. ivpSolverBenchmark

The idea around this project is first to investigate **Runge-Kutta** parallelisation across the degree of freedom of a **n-ODE system** regarding equations inter-dependancy. Our base use-case use a n*n lattice where each node are only connected to their closest neighbourgh. This induce a sparse patern in the system of equations.

Then compare those RK4 methods with a parallel-in-time algorithm such as **Parareal**. Parallelization will be tested on both CPU-Multithreading and GPU overload.

We are solving the following equation of motion:
$$ M * \ddot{X} + B * \dot{X} + K * X = F(t) $$
Were M, B and K are respectively the Mass, Damping and Stifness matrix.

# 2. Implementations

## a) Runge-Kutta 4
Implement the classic Runge-Kutta 4 algoritmh adapted to solve n-DOF second order ODE.

# 3. Make it run

1. The programm use the **Eigen/Dense** library for algebra representation and efficient BLAS: https://eigen.tuxfamily.org/index.php?title=Main_Page
2. Compile the code with the Makefile.
3. Use the Python script to run the simulation and plot the results.

# TODO
### Algorithm
- Parareal
- Dormand-Prince ?
### Optimization
- Make a better use of Eigen/Dense capability 
  - Defining butcher tab as 4(->3) Eigen Sparses Vector ?
  - Defining M, B and K as 3 Eigen Sparses Matrix
### Benchmarking
### Interface
