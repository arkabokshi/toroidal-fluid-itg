# Welcome!

This repository contains the Toroidal Fluid ITG model used in [Bokshi PPCF 2016](https://iopscience.iop.org/article/10.1088/0741-3335/58/7/075011).

Investigating the parameter space of the fluid model has been enhanced via the inclusion of a [Gaussian Process Regression (GPR)](https://github.com/C-bowman/inference-tools) library, written by Chris Bowman.

The `automation.py` script is currently set up to investigate three parameters of interest:
- `etag` 
- `epsilonn` 
- `shear`

Explanations of the above can be found in the Bokshi 2016 paper.

## Getting started

The GPR algorithm can be run by simply entering the following:

```
python3 scripts/automation.py
```

You will be prompted to input a relative file path at which to store the routine's output data, e.g:

```
./data/GPR_Test
```

Note that the output data won't be moved to this location until after the routine has finished running.

You will also be prompted to enter the number of iterations of the GPR algorithm you wish to carry out, e.g:

```
73
```

At present, the routine will first run the simulation 27 times to build an evenly distributed knowledge of the parameter space. It will then go on to iterate for as many times as requested.

A typical iteration takes ~ 20 minutes to run, so in the above example, the total run-time would be ≈ 33 hours.

Post processing can be carried out using the `gpr_plots`, `fieldprofile`, `analysis`, and `video` scripts.

Results from the GPR algorithm are saved to a csv file called `gpo_results`.

## Advanced Setup

If you wish to change the parameters of interest, you will need to make modifications to the `ToroidalFluidITG` class in `models.py`, specifically the `set_variables` method, as well as its class properties.

Care must be taken when constructing the `x` and `y` variables which are used for the instantiation of the `GpOptimiser` class in the `automation.py` script. They must both be a `numpy.ndarray`, with shape (number of points, number of dimensions) for `x`. For example, the current `x` and `y` variables look like this:

```python
x = [
        [
            eta_g_1, epsilon_n_1, shear_1
        ],
        [
            eta_g_2, epsilon_n_2, shear_2
        ],
        ...
    ]

y = [
        objective_function_value_1,
        objective_function_value_2,
        ...
    ]

# Set the bounds for [ eta_g, epsilon_n, shear ]
# as tuples of (lower_bound, upper_bound)
bounds = [
    (eta_g_lower, eta_g_upper),
    (epsilon_n_lower, epsilon_n_upper),
    (shear_lower, shear_upper)
]

# Create an instance of GpOptimiser
GPO = GpOptimiser(x, y, bounds=bounds, acquisition=UpperConfidenceBound)
```
The `GpOptimiser` class must be instantiated with a minimum of two values for `x` and `y`.
