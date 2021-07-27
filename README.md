# InternCohort2021_VASP_Plotter

To run these files:

1. Install Julia and follow all the installation instructions (if any)
2. In the terminal, navigate to the directory where you want to place these files
3. `git clone https://github.com/oashour/InternCohort2021_VASP_Plotter.git`
4. `cd` into the directory
5. `julia --project`
6. type `]` at the julia prompt. It should change to `pkg`.
7. type `instantiate`. Wait until the packages install.
8. Exit Julia using `exit()`.
9. Place a `vasprun.xml` file from whatever calculation you want (must be semiconductor with small sigma) in the folder.
10. Run `julia --project plotter.jl`
