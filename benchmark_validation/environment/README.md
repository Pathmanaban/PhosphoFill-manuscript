# Frozen benchmark environment

The completed benchmark was run in Ubuntu under WSL using:

- Python 3.12.3
- Interpreter: `/home/paddy/venvs/ptm/bin/python`
- NumPy 1.26.4
- Biopython 1.87
- pandas 2.1.4
- SciPy 1.12.0
- Matplotlib 3.8.4
- OpenMM 8.5.0

The pinned Python distributions are listed in `requirements.txt`. OpenMM may
require platform-specific installation instructions when CUDA, OpenCL, or HIP
acceleration is desired. The frozen run used the script default
`--openmm-platform auto`.
