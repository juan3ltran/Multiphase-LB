# Multiphase Lattice Boltzmann Simulation with OpenMP

This project implements a multiphase Lattice Boltzmann (LB) simulation, parallelized using OpenMP for enhanced performance. The simulation generates visual results, which can be rendered using Gnuplot.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

## Introduction

The multiphase Lattice Boltzmann simulation provides a numerical approach for simulating fluid dynamics, especially for complex multiphase flows. By utilizing OpenMP for parallelization, the computational efficiency is significantly improved.

## Features

- Simulates multiphase fluid dynamics
- Parallelized using OpenMP for enhanced performance
- Generates visual outputs compatible with Gnuplot

## Requirements

- C++ Compiler (g++)
- OpenMP library
- Gnuplot

## Installation

To compile and run the simulation, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/juan3ltran/Multiphase-LB.git
    cd Multiphase-LB
    ```

2. Compile the simulation using the provided Makefile:
    ```bash
    make
    ```

## Usage

To run the simulation and generate the output for Gnuplot, use the following commands:

1. Execute the compiled program:
    ```bash
    ./Multiphase | gnuplot
    ```

This will run the simulation and pipe the output directly to Gnuplot for rendering the results.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Contact

For any questions or suggestions, feel free to reach out:

- **Juan José** - jubeltranr@unal.edu.co / beltranjuanjo.r@gmail.com