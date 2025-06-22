# SMTDock
SMTDock is a novel approach to the protein-ligand docking problem based on [AMDock](https://github.com/Valdes-Tresanco-MS/AMDock), which employs an SMT-solver to predict the binding pose of a ligand.

## Setup

A makefile is provided to setup the system through the following commands:

1. Initial environment setup:
    ```bash
    make setup     # create the virtual environment and install dependencies
    ```

2. Activate the python environment:
    ```bash
    . .venv/bin/activate
    ```

## Usage

SMTDock can be run using the command:
```bash
python3 main.py protein ligand
```

where a related folder must already exist.

To add a new protein/ligand pair, follow the execution of AMDock using the information provided in the official [tutorial](https://github.com/Valdes-Tresanco-MS/AMDock-win/wiki).

### Example

As an example, to find the binding pose of the ligand Adenosine with the protein A2A, run the command:
```bash
python3 main.py a2a adn
```
The output of the program, including information on the new pose (in pdb format) and the related energy and runtime, can be found in the subfolder `smtdock`.

## Machine environment

The program has been tested on an Apple device; for different architectures, the setup steps may include some differences.
