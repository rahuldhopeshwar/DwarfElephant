[Mesh]
  file = meshs/simple_3layer_model.e
[]

[Variables]
active = 'temperature'
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
active = 'RBConduction'
  [./RBConduction]
    type = RBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    max_x = 0.2 #0.5
    min_x = 0.4 #0.3
    max_y = 0.2 #0.04
    min_y = 0.4 #0.02
    max_z = 0.2 #0.08
    min_z = 0.4 #0.06
   ID_first_block = 0
  [../]
[]

[BCs]
active = 'RBtop RBbottom'
  [./RBtop]
    type = RBDirichletBC
    variable = temperature
    boundary = 1
    #boundary = 2
    value = 10.00
    initial_rb_userobject = initializeRBSystem
    cache_boundaries = cacheBoundaries
    mesh_modified = false
    ID_first_block = 0
  [../]
  [./RBbottom]
    type = RBDirichletBC
    variable = temperature
    boundary = 2
    #boundary = 1
    value = 31.00
    initial_rb_userobject = initializeRBSystem
    cache_boundaries = cacheBoundaries
    mesh_modified = false
    ID_first_block = 0
  [../]
[]

[Problem]
  type = DwarfElephantRBProblem
[]

[Executioner]
  type = DwarfElephantRBSteady

  solve_type = 'Newton'
  l_tol = 1.e-8
  nl_rel_tol = 1.e-8
[]

[Functions]
active = 'cacheBoundaries'
  [./cacheBoundaries]
    type = CacheBoundaries
  [../]
[]

[UserObjects]
active = 'initializeRBSystem performRBSystem'

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystem
    parameters_filename = inputfiles/RB/ParallelLayers/TestModels/largeRBTest.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    cache_boundaries = cacheBoundaries
    offline_stage = true
    execute_on = 'initial'
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineOnlineStage

    exodus_file_name = largeRBTest

    offline_stage = true
    online_stage = true
    store_basis_functions = true

    mu_bar = 1
    online_mu = '1.05 2.5 1.05'

    execute_on = 'timestep_end'
    initial_rb_userobject = initializeRBSystem
    cache_boundaries = cacheBoundaries
  [../]
[]

[Outputs]
  print_perf_log = false
  exodus = false
  execute_on = 'timestep_end'

  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]


# ======================= Parameters for the RB method ========================

# Maximum number of basis functions that will be generated in the Offline-stage
Nmax = 20

# Name of the parameters
# Please name them mu_0, mu_1, ..., mu_n for the re-usability
parameter_names = 'mu_0 mu_1 mu_2'

# Define the minimum and maximum value of the Theta object
mu_0 = '0.01000 10.15000'
mu_1 = '0.01000 12.80000'
mu_2 = '0.01000 10.15000'

# Define the number of training sets for the Greedy-algorithm
n_training_samples = 100

# Optionally:
# Determine whether the training points are generated randomly or deterministically
deterministic_training = false

# Determine whether relative or absolute error bounds are used in the Greedy-algorithm
use_relative_bound_in_greedy = false

rel_training_tolerance = 1.e-4
#quiet_mode =  false

#normalize_rb_bound_in_greedy = true
