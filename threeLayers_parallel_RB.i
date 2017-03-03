[Mesh]
  file = parallel.msh
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

#[AuxVariables]
#  [./loadVectorF0]
#   order = FIRST
#    family = LAGRANGE
#  [../]

#  [./stiffnessMatrixA0]
#    order = FIRST
#    family = LAGRANGE
#  [../]

#  [./stiffnessMatrixA1]
#    order = FIRST
#    family = LAGRANGE
#  [../]

#  [./stiffnessMatrixA2]
#    order = FIRST
#    family = LAGRANGE
# [../]
#[]

[Kernels]
  [./RBConduction]
    type = Diffusion
    variable = temperature
 #   diag_save_in = stiffnessMatrixA0
 #   save_in = loadVectorF0
   # block = 7
  [../]

 # [./RBSandstone]
 #   type = RBDiffusion
 #   variable = temperature
 #   diag_save_in = stiffnessMatrixA1
 #   save_in = loadVectorF0
 #   block = 8
 # [../]

 # [./RBShaleBottom]
 #   type = RBDiffusion
 #   variable = temperature
 #   diag_save_in = stiffnessMatrixA2
 #   save_in = loadVectorF0
 #   block = 9
 # [../]
 []

#[Materials]
#  [./shale_top]
#    type = Shale
#    block = 7
#  [../]

#  [./sandstone]
#    type = SandStone
#    block = 8
#  [../]

#  [./shale_bottom]
#    type = Shale
#    block = 9
#  [../]
#[]

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
  [./prepareData_block7]
    type = DwarfElephantRBSystem
    system = nl0

    parameters_filename = 'threeLayers_parallel_RB.i'

    offline_stage = true
    online_stage = true
    store_basis_functions = true

    online_N = 20
    online_mu = '1.05 2.5 1.5'

    file_name = large_model

    #block = 7
  [../]
[]

[Outputs]
  exodus = true
  xda = true
  xdr = true
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

