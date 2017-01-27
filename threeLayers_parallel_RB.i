[Mesh]
  file = parallel.msh
#  type = GeneratedMesh
#  dim = 3
#  nx = 10
#  ny = 4
#  nz = 2
#  xmin = 0.0
#  xmax = 3000
#  ymin = 0.0
#  ymax = 700
#  zmin = 0.0
#  zmax = 1000
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./RBStructures]
    type = RBKernel
    variable = temperature
    muTest = 20
  [../]
[]

[Materials]
  [./shale1]
    type = Shale
    block = 7
#    block = 0
  [../]

  [./sandstone]
    type = SandStone
    block = 8
  [../]

  [./shale2]
    type = Shale
    block = 9
  [../]
[]

[BCs]
  [./bottom_shale1]
    type = DirichletBC
    variable = temperature
    boundary = 'bottom'
    value = 31
  [../]
  [./top_shale1]
    type = DirichletBC
    variable = temperature
    boundary = 'top'
    value = 10
  [../]
[]

#[UserObjects]
#[]

[RBSimpleConstruction]
[]


[Executioner]
  type = Steady
  solve_type = 'PJFNK'
[]

[Outputs]
  #exodus = true

  [./Console]
    type = Console
    #execute_on = 'final'
  [../]

  [./RBOutput]
    type = RBOutput

    execute_on = 'timestep_end'
    execute_system_information_on = 'timestep_end'
    execute_postprocessors_on = 'timestep_end'

    parameters_filename = '/home/bl1/projects/moose/Conduction/comparison_different_layer_types/threeLayers_parallel_RB.i'

    offline_stage = true
    online_stage = true
    store_basis_functions = true

    online_N = 20
    online_mu0= 1.0
  [../]
[]


# ====================== Parameters for the RB approximation ======================

# Maximum number of basis functions that will be generated in the Offline-stage
Nmax = 20

# Name of the parameters
# Please name them mu_0, mu_1, ..., mu_n for the re-usability
parameter_names = 'mu_0'

# Define the minimum and maximum value of the Theta object
mu_0 = '1.1 3'

# Define the number of training sets for the Greedy-algorithm
n_training_samples = 100

# Optionally:
# Determine whether the training points are generated randomly or deterministically
deterministic_training = true

# Determine whether relative or absolute error bounds are used in the Greedy-algorithm
use_relative_bound_in_greedy = false

