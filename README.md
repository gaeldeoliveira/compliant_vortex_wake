# dogoro_2d_body

Solution of the Euler equations in vorticity form for actuation surface flows - version for planar flow with bodies

# What is this repository for ?

Collaboration between Gael de Oliveira, Vinit Dighe and Francesco Avallone.

# Where did the root code come from ?

From dogoro_axi/09_new_solvers/kirikou_HS_coupled_solver, which itselft combined the kirikou code for actuator disks (single actuator, variable loading distribution version) with the Acti_HS panel code for collections of lifting bodies.

# What is being done on the code right now ?

We are implementing shape morphing and mesh generation capabilities, to couple the Euler solver with:
-> Scripts for parametric studies
-> Optimization algorithms
-> RANS solvers like Fluent and, eventually, OpenFoam and REFRESCO

# What should and can be done, but is moving slowly ?