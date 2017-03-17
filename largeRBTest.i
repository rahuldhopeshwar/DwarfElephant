[Mesh]
  file = parallel_test.msh
[]

[Variables]
  [./temperature]
  [../]
[]

[AuxVariables]
active = 'A0'
  [./A0]
  [../]

  [./A1]
  [../]

  [./A2]
  [../]
[]

[Kernels]
  [./RBConduction_block0]
    type = RBDiffusion
    variable = temperature
    diag_save_in = A0
    block = 0
  [../]

  [./RBConduction_block1]
    type = RBDiffusion
    variable = temperature
   # diag_save_in = A1
    block = 1
  [../]

  [./RBConduction_block2]
    type = RBDiffusion
    variable = temperature
   # diag_save_in = A2
    block = 2
  [../]
 []

[Materials]
active = ''

  [./shale]
    type = Shale
    block = 0
  [../]

  [./sandstone]
    type = SandStone
    block = 1
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

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_rest'
  petsc_options_value = 'hypre  boomeramg   101'
[]

[UserObjects]
active = 'performRBSystem'

 [./prepareData_block0]
   type = DwarfElephantPrepareRBSystem
   block = 0
 [../]

 [./performRBSystem]
    type = DwarfElephantRBSystem

    parameters_filename = 'largeRBTest.i'

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
  print_perf_log = true
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

quiet_mode =  false

