# Step 2: Magma inflation with evolution of porosity

```{include} step02_inflation-synopsis.md
```

## Simulation parameters

We extend the simulation in Step 1 by including evolution of the porosity, which depends on the time derivative of the pressure and trace strain.
We also compute the deformation relative to a uniform reference compressive pressure of 5 MPa to illustrate how to use a reference state with poroelasticity.
We use the same initial conditions and boundary conditions as in Step 1.

Because the evolution of porosity depends on the time derivative of the solution subfields, we need to include the time derivatives in the solution field.
As a result, we have 6 subfields in our solution field.

```{code-block} cfg
---
caption: Solution subfields for Step 2.
---
# Poroelasticity with porosity state variable requires solution with time derivaties 
solution = pylith.problems.SolnDispPresTracStrainVelPdotTdot

# Set basis order for all solution subfields
[pylithapp.problem.solution.subfields]
displacement.basis_order = 2
pressure.basis_order = 1
trace_strain.basis_order = 1
velocity.basis_order = 2
pressure_t.basis_order = 1
trace_strain_t.basis_order = 1
```

```{code-block} cfg
---
caption: Material parameters for poroelasticity with state variables and reference state for Step 2.
---
[pylithapp.problem.materials.crust]
use_state_variables = True

db_auxiliary_field.values = [
    solid_density, fluid_density, fluid_viscosity, porosity, shear_modulus, drained_bulk_modulus, biot_coefficient, fluid_bulk_modulus, solid_bulk_modulus, isotropic_permeability,
    reference_stress_xx, reference_stress_yy, reference_stress_zz, reference_stress_xy,
    reference_strain_xx, reference_strain_yy, reference_strain_zz, reference_strain_xy
    ]
db_auxiliary_field.data   = [
    2500*kg/m**3, 1000*kg/m**3, 0.001*Pa*s, 0.01, 6.0*GPa, 10.0*GPa, 1.0, 2.0*GPa, 20.0*GPa, 1e-15*m**2,
    -5.0*MPa, -5.0*MPa, -5.0*MPa, 0.0*MPa,
    0.0, 0.0, 0.0, 0.0
    ]

auxiliary_subfields.porosity.basis_order = 1

[pylithapp.problem.materials.crust.bulk_rheology]
use_reference_state = True


[pylithapp.problem.materials.intrusion]
use_state_variables = True

db_auxiliary_field.values = [
    solid_density, fluid_density, fluid_viscosity, porosity, shear_modulus, drained_bulk_modulus, biot_coefficient, fluid_bulk_modulus, solid_bulk_modulus, isotropic_permeability,
    reference_stress_xx, reference_stress_yy, reference_stress_zz, reference_stress_xy,
    reference_strain_xx, reference_strain_yy, reference_strain_zz, reference_strain_xy
    ]
db_auxiliary_field.data   = [
    2500*kg/m**3,  1000*kg/m**3, 0.001*Pa*s, 0.1, 6.0*GPa, 10.0*GPa, 0.8, 2.0*GPa, 20.0*GPa, 1e-13*m**2,
    -5.0*MPa, -5.0*MPa, -5.0*MPa, 0.0*Pa,
    0.0, 0.0, 0.0, 0.0
    ]

auxiliary_subfields.porosity.basis_order = 1

[pylithapp.problem.materials.intrusion.bulk_rheology]
use_reference_state = True
```

## Running the simulation

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_inflation.cfg

# The output should look something like the following.

# -- many lines ommitted --

```

In contrast with Step 1, the evolution of porosity in Step 2 results in nonlinear governing equations.
This requires multiple iterations of the linear solve for the nonlinear solver to converge.

## Visualizing the results

In {numref}`fig:example:magma:2d:step02:solution` we use ParaView to visualize the evolution of the y displacement component using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/magma-2d` directory.

```{code-block} console
---
caption: Open ParaView using the command line.
---
$ PATH_TO_PARAVIEW/paraview

# For macOS, it will be something like
$ /Applications/ParaView-5.10.1.app/Contents/MacOS/paraview
```

Next we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.
For Step 2 we need to change the simulation.

:::{figure-md} fig:example:magma:2d:step01:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2 at t=100 yr. The colors indicate the fluid pressure, and the deformation is exaggerated by a factor of 1000." width="75%"/>

Solution for Step 2 at t=100 yr.
The colors of the shaded surface indicate the fluid pressure, and the deformation is exaggerated by a factor of 1000.
:::
