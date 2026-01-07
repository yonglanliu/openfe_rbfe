# OpenFE RBFE Dashboard

An intuitive, web-based dashboard for performing **Relative Binding Free Energy (RBFE)** and **alchemical Free Energy Perturbation (FEP)** calculations using the [OpenFE framework](https://github.com/IIIS-Li-Group/OpenFE).

This project bridges the gap between powerful alchemical free energy methods and practical usability by providing an end-to-end graphical interface for system setup, execution, and analysis.

---

## 📖 Background & Concepts

If you are new to alchemical transformations, start with this introductory guide:

👉 **Read the Blog Post:**  
[Introduction to Alchemical Free Energy Calculation: FEP and TI for RBFE](https://yonglanliu.github.io/2025/12/19/RBFE.html)

The post explains:
- What RBFE and FEP are
- When and why to use alchemical methods
- Conceptual intuition behind alchemical transformations  
—all without diving deep into heavy mathematics.

---

## 🚀 Key Features

This dashboard streamlines the standard OpenFE RBFE workflow into **five user-friendly steps**:

1. **Molecular Visualization**  
   - Interactively view protein targets and ligand structures
   - Inspect binding poses and ligand alignments directly in the browser

2. **Network Planning**  
   - Automatically generate ligand transformation networks
   - Visualize and inspect alchemical edges
   - Identify disconnected components or problematic mappings

3. **System Setup**  
   - Configure solvation models, force fields, and simulation boxes
   - Define alchemical protocols and λ schedules
   - Validate system integrity before job submission

4. **Job Execution**  
   - Prepare OpenFE-compatible simulation scripts
   - Launch and monitor simulations from a unified interface
   - Designed to integrate with local, HPC, or cloud-based execution

5. **Results Analysis**  
   - Parse raw OpenFE output data
   - Visualize convergence and overlap
   - Compute final **ΔΔG** values with uncertainty estimates

---

## 🧠 Design Philosophy

- **Reproducibility first** – all steps are transparent and configurable  
- **Minimal boilerplate** – reduce repetitive scripting  
- **Beginner-friendly** – sensible defaults with advanced overrides  
- **Power-user capable** – direct access to OpenFE primitives when needed  

---

## 🛠️ Technology Stack

- **Backend:** Python, OpenFE
- **Frontend:** Web-based UI (framework-dependent)
- **Visualization:** Molecular viewers for proteins and ligands
- **Data Processing:** OpenFE analysis utilities

> ⚠️ Exact frontend/backend frameworks may evolve as the project matures.

---

## 📦 Installation

### Prerequisites

- Python ≥ 3.9
- Conda or Mamba (recommended)
- OpenFE and its dependencies
- A working OpenMM-compatible environment

### Clone the Repository

```bash
git clone https://github.com/your-username/openfe-rbfe-dashboard.git
cd openfe-rbfe-dashboard
```
### Install
```bash
conda env create -f environment.yml
conda activate openfe-dashboard
```

### Run
```bash
stremlit run Home.py
```
