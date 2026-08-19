<p align="center">
  <img src="assets/flow.svg" alt="Potential flow past a cylinder, resolved daily" width="100%">
</p>

<h1 align="center">Tyler Jones</h1>

<p align="center">
  <b>CFD Engineer @ Cadence</b> &nbsp;·&nbsp; <b>M.S. student @ UW–Madison</b><br>
  Finite volume methods, solver development, and the occasional 18M-cell mesh that refuses to converge.
</p>

<p align="center">
  <a href="https://www.linkedin.com/in/YOUR-HANDLE"><img src="https://img.shields.io/badge/LinkedIn-0A66C2?style=flat-square&logo=linkedin&logoColor=white" alt="LinkedIn"></a>
  <a href="mailto:YOUR-EMAIL"><img src="https://img.shields.io/badge/Email-D14836?style=flat-square&logo=gmail&logoColor=white" alt="Email"></a>
  <img src="https://komarev.com/ghpvc/?username=tjjones6&style=flat-square&color=blue" alt="Profile views">
</p>

---

### About

I build and validate computational fluid dynamics solvers — both commercially and from scratch. Day to day I run high-lift validation cases (CRM-HL / AIAA High Lift Prediction Workshop) in Cradle scFLOW and Fidelity CFD. Evenings and weekends I write my own solvers in C++ and Python to understand what the commercial codes are actually doing under the hood.

- 🔬 **Currently building:** `cavity2d` — a from-scratch C++ incompressible Navier–Stokes solver (SIMPLE + Rhie–Chow on a collocated grid) with geometric multigrid acceleration of the pressure-correction equation
- 📚 **Studying:** numerical methods for PDEs, Krylov solvers, and preconditioning
- 🧩 **Interested in:** high-order schemes, turbulence modeling (SA, k-ω SST), HPC, and making solvers that are actually readable
- 💬 **Ask me about:** OpenFOAM mesh conversion pain, Rhie–Chow interpolation, or why your residuals plateaued at 1e-3

---

### Selected Work

| Project | Description | Stack |
| :--- | :--- | :--- |
| **[cavity2d](https://github.com/tjjones6)** | Incompressible NS solver: collocated FVM, SIMPLE, Rhie–Chow, multigrid pressure correction. Verified via manufactured solutions, validated against Ghia et al. (1982). | `C++17` `CMake` `VTK` |
| **[naca-viscous](https://github.com/tjjones6)** | Full viscous–inviscid airfoil analysis: panel method + integral boundary layer + coupling. XFOIL, but mine. | `Python` `NumPy` |
| **[sa-channel](https://github.com/tjjones6)** | Spalart–Allmaras turbulent channel flow solver written from the model equations up. | `Python` |
| **[hlpw-tools](https://github.com/tjjones6)** | Boundary-condition calculators and log parsers for High Lift Prediction Workshop cases. | `Python` `Tkinter` |

<details>
<summary><b>📖 A bit more on the multigrid work</b></summary>

<br>

The pressure-correction (Poisson) solve dominates cost in a SIMPLE loop — a plain Gauss–Seidel smoother converges at a rate that degrades as `O(h²)` with mesh refinement, because low-frequency error components are essentially invisible to a local relaxation scheme.

Geometric multigrid fixes this by restricting the residual onto coarser grids where those smooth modes *are* high-frequency, solving there, and prolongating the correction back. Done right, the iteration count becomes mesh-independent.

The interesting engineering question is where the V-cycle stops paying for itself on a collocated grid where Rhie–Chow momentum interpolation couples the pressure and velocity fields in a way that isn't quite the clean Laplacian the textbooks assume.

</details>

<details>
<summary><b>🛠️ Tooling & environment</b></summary>

<br>

**Solvers:** OpenFOAM v2312 · Cradle scFLOW · Fidelity CFD
**Languages:** C++ · Python · MATLAB · Bash
**Build & VC:** CMake · Git · wmake
**Meshing:** Pointwise · Gmsh · blockMesh (usually via a Python generator)
**Post:** ParaView · VTK · Matplotlib
**Environment:** WSL2 / Ubuntu, because the Windows laptop was not negotiable

</details>

---

### Stats

<p align="center">
  <img height="165" src="https://github-readme-stats.vercel.app/api?username=tjjones6&show_icons=true&hide_border=true&theme=default&include_all_commits=true" alt="GitHub stats">
  <img height="165" src="https://github-readme-stats.vercel.app/api/top-langs/?username=tjjones6&layout=compact&hide_border=true&theme=default" alt="Top languages">
</p>

---

<p align="center"><i>"Essentially, all models are wrong, but some are useful." — George Box</i></p>
