[Mesh]
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

  block_id = '0 1'
  block_name = 'shale sandstone'
[]

[MeshModifiers]
  [./subdomains]
    type = AssignElementSubdomainID
  #element_id = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 39 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 646 56 66 67 68 69 70 71 72 73 74 75 76 77 78 79'
  subdomain_ids = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'
  [../]
[]

[Variables]
active = 'temperature'
  [./temperature]
  [../]
[]

[Kernels]
active = 'RBConduction_block0 RBConduction_block1'
  [./RBConduction_block0]
    type = RBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    block = 0
    subdomain = 0
  [../]

  [./RBConduction_block1]
    type = RBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    block = 1
    subdomain = 1
  [../]
 []

[Materials]
active = ''

  [./shale_top]
    type = Shale
    block = 2
  [../]

  [./sandstone]
    type = SandStone
    block = 1
  [../]

  [./shale_bottom]
    type = Shale
    block = 0
  [../]
[]

[BCs]
active = 'bottom top'
  [./bottom]
    type = RBDirichletBC
    variable = 'temperature'
    boundary = 'bottom'
    value = 31
    initial_rb_userobject = initializeRBSystem
    cache_stiffness_matrix = cacheStiffnessMatrix
  [../]

  [./top]
    type = RBDirichletBC
    variable = 'temperature'
    boundary = 'top'
    value = 10
    initial_rb_userobject = initializeRBSystem
    cache_stiffness_matrix = cacheStiffnessMatrix
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_rest'
  petsc_options_value = 'hypre  boomeramg   101'
[]

[Functions]
  [./cacheStiffnessMatrix]
    type = CacheStiffnessMatrix
  [../]
[]

[UserObjects]
active = 'initializeRBSystem performRBSystem'

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystem
    parameters_filename = smallRBTest.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    offline_stage = true
    online_stage = true
    store_basis_functions = true
    execute_on = initial
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineStage

    parameters_filename = smallRBTest.i
    residual_name = Re_non_time

    offline_stage = true
    online_stage = true
    store_basis_functions = true

    mu_bar = 1
    online_N = 2
    online_mu = '1.05 2.5' # 1.15'

    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true

    execute_on = timestep_end
    initial_rb_userobject = initializeRBSystem
    cache_stiffness_matrix = cacheStiffnessMatrix
  [../]
[]

[Outputs]
  exodus = false
  execute_on = 'timestep_end'
#  print_perf_log = true
[]

# ====================== Parameters for the RB approximation ======================

# Maximum number of basis functions that will be generated in the Offline-stage
Nmax = 20

# Name of the parameters
# Please name them mu_0, mu_1, ..., mu_n for the re-usability
parameter_names = 'mu_0 mu_1' # mu_2'

# Define the minimum and maximum value of the Theta object
mu_0 = '1.99 2.01'
mu_1 = '2.2 2.8'
#mu_2 = '0.95 1.15'

# Define the number of training sets for the Greedy-algorithm
n_training_samples = 100

# Optionally:
# Determine whether the training points are generated randomly or deterministically
deterministic_training = true

# Determine whether relative or absolute error bounds are used in the Greedy-algorithm
use_relative_bound_in_greedy = true

rel_training_tolerance = 1e-4
quiet_mode =  false
#evaluate_RB_error_bound = false

#normalize_rb_bound_in_greedy = true
