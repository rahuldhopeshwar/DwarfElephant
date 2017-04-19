[Mesh]
  file = parallel.msh
  block_id = '7 8 9'
  block_name = 'shale_top sandstone shale_bottom'
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./RBConduction]
    type = Diffusion
    variable = temperature
  [../]
 []

[Materials]
active = ''
  [./shale_top]
    type = Shale
    block = 7
  [../]

  [./sandstone]
    type = SandStone
    block = 8
  [../]

  [./shale_bottom]
    type = Shale
    block = 9
  [../]
[]

[BCs]
  [./bottom]
    type = DirichletBC
    variable = temperature
    boundary = 'bottom'
    value = 31
  [../]
  [./top]
    type = DirichletBC
    variable = temperature
    boundary = 'top'
    value = 10
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
[]

[UserObjects]
active = 'performRBSystem'

  [./prepareData_block7]
    type = DwarfElephantPrepareRBSystem
    block = 7
  [../]

  [./prepareData_block8]
    type = DwarfElephantPrepareRBSystem
    block = 8
  [../]

  [./prepareData_block9]
    type = DwarfElephantPrepareRBSystem
    block = 9
  [../]

  [./performRBSystem]
    type = DwarfElephantRBSystem

    parameters_filename = 'threeLayer_parallel_RB.i'

    offline_stage = true
    online_stage = true
    store_basis_functions = true

    online_N = 20
    online_mu = '1.05 2.5 1.5'
  [../]
[]

[Outputs]
  exodus = true
  execute_on = 'timestep_end'
  #print_perf_log = true
[]

# ====================== Parameters for the RB approximation ======================

# Maximum number of basis functions that will be generated in the Offline-stage
Nmax = 20

# Name of the parameters
# Please name them mu_0, mu_1, ..., mu_n for the re-usability
parameter_names = 'mu_0 mu_1 mu_2'

# Define the minimum and maximum value of the Theta object
mu_0 = '0.95 1.15'
mu_1 = '2.2 2.8'
mu_2 = '0.95 1.15'

# Define the number of training sets for the Greedy-algorithm
n_training_samples = 10

# Optionally:
# Determine whether the training points are generated randomly or deterministically
deterministic_training = false

# Determine whether relative or absolute error bounds are used in the Greedy-algorithm
use_relative_bound_in_greedy = false

quiet_mode = false
