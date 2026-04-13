# pcebatt (PCE-based Battery Cell Simulation)

This repository contains the pcebatt Python package for battery cell simulation.

## Prerequisites
Before installing, ensure you have Python 3.10 or higher. It is highly recommended to use a virtual environment (Conda or venv) to avoid dependency conflicts.

## Installation
### 1. Clone the Repository
Clone this repository to your local machine and navigate into the folder:
```Bash
git clone https://github.com/AregM02/PCE_Batt.git
cd PCE_Batt
```

### 2. Install the Package
Once you are in the project directory, you can install the package:

```Bash 
pip install .
```

### 3. Verification
To ensure the package is installed correctly and all dependencies are resolved, run the following line. This confirms Python is using the globally installed version in site-packages rather than the local source folder:

```Bash
python -c "import pcebatt; print(f'Installed at: {pcebatt.__file__}')"
```

## Running Tests
To verify proper functionality, you can use the included tests. To do so, follow the steps below.

<ins>Install test dependencies</ins>: Ensure you are inside the PCE_Batt directory and run:
```Bash
pip install ".[test]"
```

<ins>Run the validation test</ins>:
Navigate away from the repository folder (e.g., up one level) and point pytest to the test script:
```Bash
cd ..
pytest PCE_Batt/tests/test_simulation.py
```
## Usage
You can find examples of how to utilize the pcebatt package within the examples/ directory. These notebooks demonstrate how to configure simulations and interpret results using the PCE framework.

<ins>Project Structure</ins>

    pcebatt/: The core source code for the battery simulation package.

    examples/: Jupyter notebooks and scripts demonstrating package usage.

<ins>Troubleshooting</ins>

    Dependency Conflicts: If you encounter issues with numpy or scipy, ensure your virtual environment is clean and that you are using a compatible version of Python (>=3.10).

    Missing Data: The MANIFEST.in file is configured to include all files within the pcebatt/ directory. If you have added new data files to the package and they are not being included in the build, ensure they are tracked by Git and that the MANIFEST.in rules cover them.
