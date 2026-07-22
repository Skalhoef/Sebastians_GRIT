# GRIT

### Green Functions for Interacting Thermodynamics

GRIT is a parallel Fortran research code for calculating finite-temperature
thermodynamic properties of crystalline materials from first-principles
electronic, phononic, and electron–phonon data. It sits at the end of a
Quantum ESPRESSO–Wannier90–EPW workflow and evaluates both non-interacting
and electron–phonon-interacting contributions using Green-function methods.

The code was developed for computational condensed-matter research, with an
emphasis on large Brillouin-zone datasets and execution on HPC systems.

> **Project status:** Active research software. GRIT is configured in the
> source and expects precomputed datasets that are not distributed with this
> repository. The interface and file formats may evolve.

## What GRIT computes

- Electronic state counting and particle-number consistency at $T=0$
- Temperature-dependent chemical potential at fixed particle number
- Specific heat per unit cell, with optional electron–phonon corrections
- Electronic and phononic spectral contributions from Green functions
- Optional refinement around the Fermi surface and the Brillouin-zone centre
- Conversion of large text datasets to binary files for faster subsequent I/O

The source contains material presets for Ag, Au, C (diamond), Cs, Cu, K, Pb,
Rb, and Si. These presets define quantities such as the Fermi energy, band
count, sampling grids, broadening, and temperature range; they are examples
rather than bundled reference calculations.

## Computational workflow

```mermaid
flowchart LR
    QE[Quantum ESPRESSO] --> W90[Wannier90]
    W90 --> EPW[EPW / WannierTools]
    EPW --> DATA[Fine-grid energies, phonons, and couplings]
    DATA --> GRIT[GRIT]
    GRIT --> OUT[Chemical potential and specific heat vs temperature]
```

GRIT consumes the following quantities:

1. The Fermi energy from a Quantum ESPRESSO NSCF calculation
2. Electronic energies on a fine $k$-point grid from EPW or WannierTools
3. Phonon energies on a fine $q$-point grid from Quantum ESPRESSO or EPW
4. Electron–phonon matrix elements on fine grids from EPW when interacting
   calculations are enabled

The upstream EPW version used during development includes additional output
options for exporting fine-grid data. It is available in the related
[Quantum ESPRESSO and EPW repository](https://github.com/Skalhoef/Sebastians_Personalized_QE_n_EPW_Codes).
The relevant EPW flags are:

```fortran
print_fine_Fermi = .true.
prtgkk_sebbe     = .true.
```

## Parallel design

GRIT is MPI-parallel and is built with OpenMP support:

- MPI distributes fine-grid $k$-point slices between processes.
- Each process reads only its local slice of the largest coupling arrays.
- OpenMP runtime support and scheduler settings are present as scaffolding for
  thread-level parallelism; the current computational distribution uses MPI.
- Stream-format binary caches reduce the cost of repeatedly loading large
  text exports.

The repository includes build configurations for a workstation and for the
Dardel supercomputer at PDC.

## Requirements

- A Fortran compiler with Fortran 2008 support
- An MPI implementation that provides the `mpif90` wrapper
- OpenMP support
- GNU Make
- Precomputed Quantum ESPRESSO/EPW data in the layout expected by
  `set_filepaths()` in `modmain.f90`

The workstation Makefile is written for GNU-style compiler flags. If your MPI
wrapper uses another compiler, adjust `FC` and `FFLAGS` in `Makefile_PC`.

## Build

```bash
make -f Makefile_PC
```

Other useful targets are:

```bash
make -f Makefile_PC clean
make -f Makefile_PC rebuild
make -f Makefile_PC help
```

This produces the executable `grit.x`.

For Dardel, load the site-specific compiler environment and build with:

```bash
make -f Makefile_Dardel
```

## Configure a calculation

GRIT currently uses compile-time configuration. Before building, edit the
input-parameter section near the top of `modmain.f90`:

1. Select one `structure_name` or add a material in
   `set_structure_parameters()`.
2. Choose the calculation modes: `interacting`, `Fermi_centered_k`, and
   `Gamma_centered_q`.
3. Enable the required tasks, such as `calc_mu_vs_T` or `calc_c_V_per_UC`.
4. Review the temperature range, energy-grid resolution, broadenings, and
   numerical thresholds.
5. Update `set_filepaths()` if your data directory differs from the expected
   project layout.

The default path convention is:

```text
../B_QE_n_EPW_Files/Step3_EPW/
└── <material>/
    └── Dispersions_<Nk>_<Nq>_<Nq_SE>_<fsthick>_<Ni_uniform>/
        ├── epsilon_data.txt
        ├── omega_data.txt
        ├── SE_epsilon_data.txt
        ├── SE_epsilon_kplusq_data.txt
        ├── SE_g_data.txt
        ├── SE_omega_data.txt
        ├── uniform_epsilon_data.txt
        ├── uniform_g_data.txt
        ├── uniform_k_data.txt
        └── uniform_omega_data.txt
```

Only the files required by the selected modes and tasks need to be present.

## Run

Run locally with any suitable number of MPI processes:

```bash
mpirun -np 4 ./grit.x
```

On Dardel, `job_Fortran.sh` provides a SLURM template. Review the allocation,
partition, resource request, and module versions before submitting it with
`sbatch`.

## Outputs

Depending on the enabled tasks, GRIT writes plain-text results alongside the
input data, including:

- `GRIT_chemical_potential_vs_T_*.txt`
- `GRIT_c_V_per_uc_vs_T_*.txt`
- `GRIT_temperature_values_*.txt`
- Energy-resolved diagnostic integrands when debug mode is enabled

When `create_binary` is enabled for an interacting calculation, GRIT also
creates binary counterparts of the largest coupling and $\epsilon_{k+q}$
datasets. These caches are compiler/platform-dependent unformatted data and
should be regenerated when moving between incompatible systems.

## Repository structure

| Path | Purpose |
| --- | --- |
| `grit.f90` | Program entry point and task orchestration |
| `modmain.f90` | Configuration, physical constants, shared state, and utilities |
| `utils.f90` | Numerical kernels, Green functions, MPI distribution, and data I/O |
| `Makefile_PC` | Workstation build using `mpif90` |
| `Makefile_Dardel` | PDC Dardel build using the Cray `ftn` wrapper |
| `job_Fortran.sh` | Example hybrid MPI/OpenMP SLURM job |

## Current limitations

- Calculation settings and material definitions are compiled into the binary.
- Input datasets and automated regression tests are not yet included.
- The Sommerfeld-coefficient task (`calc_gamma`) is experimental.
- Binary caches are intended for local performance, not portable data exchange.

## Reproducibility and contributions

For reproducible work, record the GRIT commit, compiler and MPI versions,
material preset, enabled modes, numerical parameters, and provenance of every
upstream dataset. Issues and focused pull requests are welcome; please include
a minimal description of the configuration and data dimensions involved.

## Author

Developed by [Sebastian Skalhoef](https://github.com/Skalhoef) as part of a
computational condensed-matter physics workflow.
