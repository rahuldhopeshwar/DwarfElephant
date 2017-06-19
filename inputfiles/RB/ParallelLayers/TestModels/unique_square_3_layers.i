[Mesh]
  file = meshs/unit_3layer.e
  #file = meshs/simple_3_model_unit_volume_global_mesh_factor_0-1.msh
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
#active = 'Conduction'
#active = 'Conduction Euler'
  [./RBConduction]
    type = RBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    simulation_type = transient
    #vector_seperation_according_to_subdomains = false
  [../]

  [./Conduction]
    type = Conduction
    variable = temperature
  [../]

  [./Euler]
    type = TimeDerivative
    variable = temperature
  [../]
[]

[BCs]
active = 'RBtop RBbottom'
#active = 'top bottom'
  [./RBtop]
    type = RBDirichletBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 2 #4
    value = 10.00
    initial_rb_userobject = initializeRBSystem
    cache_boundaries = cacheBoundaries
    mesh_modified = false
    simulation_type = transient
    ID_Fq = 2
    ID_Aq = 2
  [../]
  [./RBbottom]
    type = RBNeumannBC
    variable = temperature
    boundary = 1 #2
    value = -38.96
    cache_boundaries = cacheBoundaries
    initial_rb_userobject = initializeRBSystem
    mesh_modified = false
    simulation_type = transient
    ID_Fq = 0
    ID_Aq = 0
  [../]

  [./top]
    type = DirichletBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 2
    value = 10.00
  [../]
  [./bottom]
    type = NeumannBC
    variable = temperature
    #boundary = 'leftbottom rightbottom'
    boundary = 1
    value = 27.275
  [../]
[]

[Materials]
active = ' '
#active = 'shale_top sandstone shale_bottom'
  [./shale_top]
    type = Shale
    block = 0
  [../]

  [./sandstone]
    type = SandStone
    block = 1
  [../]

  [./shale_bottom]
    type = Shale
    block = 2
  [../]
[]

[Problem]
  type = DwarfElephantRBProblem
[]

[Executioner]
  type = DwarfElephantRBSteady
  #type = Steady
  #type = Transient

  #num_steps = 35
  #dt = 0.125

  solve_type = 'Newton'
  l_tol = 1.e-8
  nl_rel_tol = 1.e-8
[]

[Functions]
active = 'cacheBoundaries'
#active = ''
  [./cacheBoundaries]
    type = CacheBoundaries
  [../]
[]

[UserObjects]
active = 'initializeRBSystem performRBSystem'
#active = ''

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystemTransient
    parameters_filename = inputfiles/RB/ParallelLayers/TestModels/unique_square_3_layers.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    cache_boundaries = cacheBoundaries
    offline_stage = true
    execute_on = 'initial'
    transient = true
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineOnlineStageTransient

    exodus_file_name = unique_square_3_layers

    offline_stage = true
    online_stage = true
    offline_error_bound = true
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
mu_0 = '1.01000 5.15000'
mu_1 = '1.01000 7.15000'
mu_2 = '1.01000 5.15000'

# Define the number of training sets for the Greedy-algorithm
n_training_samples = 100

# Optionally:
# Determine whether the training points are generated randomly or deterministically
deterministic_training = false

# Determine whether relative or absolute error bounds are used in the Greedy-algorithm
use_relative_bound_in_greedy = false

rel_training_tolerance = 1.e-5
#quiet_mode =  false

#normalize_rb_bound_in_greedy = true
