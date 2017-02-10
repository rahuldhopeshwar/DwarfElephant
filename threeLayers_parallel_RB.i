[Mesh]
#  file = parallel.msh
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 4
  nz = 2
  xmin = 0.0
  xmax = 3000
  ymin = 0.0
  ymax = 700
  zmin = 0.0
  zmax = 1000
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE

  [../]
[]

[AuxVariables]
  [./loadVectorF0]
    order = FIRST
    family = LAGRANGE
  [../]

  [./stiffnessMatrixA0]
    order = FIRST
    family = LAGRANGE
  [../]

  [./stiffnessMatrixA1]
    order = FIRST
    family = LAGRANGE
  [../]

  [./stiffnessMatrixA2]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./RBShaleTop]
    type = RBDiffusion
    variable = temperature
    diag_save_in = stiffnessMatrixA0
    save_in = loadVectorF0
    block = 0
  [../]

#  [./RBSandstone]
#    type = RBDiffusion
#    variable = temperature
#    diag_save_in = stiffnessMatrixA1
#    save_in = loadVectorF0
#    block = 8
#  [../]

#  [./RBShaleBottom]
#    type = RBDiffusion
#    variable = temperature
#    diag_save_in = stiffnessMatrixA2
#    save_in = loadVectorF0
#    block = 9
#  [../]
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

[Outputs]
  exodus = true
  xda = true
  xdr = true
  execute_on = 'timestep_end'

  [./RBOutput]
    type = RBOutput

    parameters_filename = '/home/bl1/projects/DwarfElephant/threeLayers_parallel_RB.i'

    offline_stage = true
    online_stage = true
    store_basis_functions = true

    online_N = 20
    online_mu0 = 1.05
  [../]

  [./F0]
    type = XDR
    show = loadVectorF0
    execute_on = 'timestep_end'
    file_base = loadVectorF0
  [../]

  [./A0]
    type = XDR
    show = stiffnessMatrixA0
    execute_on = 'timestep_end'
    file_base = stiffnessMatrixA0
  [../]

  [./A1]
    type = XDR
    show = stiffnessMatrixA1
    execute_on = 'timestep_end'
    file_base = stiffnessMatrixA1
  [../]

  [./A2]
    type = XDR
    show = stiffnessMatrixA2
    execute_on = 'timestep_end'
    file_base = stiffnessMatrixA2
  [../]
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
n_training_samples = 100

# Optionally:
# Determine whether the training points are generated randomly or deterministically
deterministic_training = false

# Determine whether relative or absolute error bounds are used in the Greedy-algorithm
use_relative_bound_in_greedy = false

