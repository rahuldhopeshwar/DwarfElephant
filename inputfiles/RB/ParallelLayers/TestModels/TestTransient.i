[Mesh]
#file = meshs/unit_2layer.e
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmin = 0.0
  xmax = 1
  ymin = 0.0
  ymax = 1
  zmin = 0.0
  zmax = 1
#  #elem_type=TET4
[]

[MeshModifiers]
active = 'subdomains'
#active = ' '
  [./subdomains]
    type = AssignElementSubdomainID
    subdomain_ids = '0 0 1 1 0 0 1 1'
  [../]
[]

[Variables]
active = 'temperature'
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
active = 'RBConduction RBTime'
#active = 'Conduction'
#active = 'Conduction Euler'
  [./RBConduction]
    type = DwarfElephantRBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    simulation_type = transient
  [../]

  [./RBTime]
    type = DwarfElephantRBTimeDerivative
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    simulation_type = transient
  [../]

  [./Conduction]
    type = DwarfElephantConduction
    variable = temperature
    lifting_function = temperature_gradient
    initial_rb_userobject = initializeRBSystem
  [../]

  [./Euler]
    type = TimeDerivative
    variable = temperature
  [../]
[]

[BCs]
active = 'RBtop RBbottom'
#active = 'RBtop'
#active = 'top bottom'
#active = ' '
  [./RBtop]
    type = DwarfElephantRBPresetBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 3 #3 (MOOSE)
    value = 0.00
    initial_rb_userobject = initializeRBSystem
    ID_Aq = 1
    simulation_type = transient
  [../]
  [./RBbottom]
    type = DwarfElephantRBNeumannBC
    variable = temperature
    boundary = 1 #2
    value = -40
    initial_rb_userobject = initializeRBSystem
    ID_Aq = 0
    simulation_type = transient
  [../]

  [./top]
    type = DirichletBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 3
    value = 0
  [../]
  [./bottom]
    type = NeumannBC
    variable = temperature
    #boundary = 'leftbottom rightbottom'
    boundary = 1
    value = 40
  [../]
  [./left]
    type = FunctionDirichletBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 5
    function = temperature_gradient
  [../]
[]

[Materials]
active = ' '
#active = 'shale sandstone'
  [./shale]
    type = DwarfElephantShale
    block = 0
  [../]

  [./sandstone]
    type = DwarfElephantSandStone
    block = 1
  [../]
[]

[Problem]
  type = DwarfElephantRBProblem
[]

[Executioner]
  type = DwarfElephantRBExecutioner
  #type = Steady
  #type = Transient

  #num_steps = 10
  #dt = 0.01

  solve_type = 'Newton'
  l_tol = 1.e-8
  nl_rel_tol = 1.e-8
[]

[Functions]
#active = 'cacheBoundaries'
#active = 'temperature_gradient'
active = ' '
  [./cacheBoundaries]
    type = CacheBoundaries
  [../]

  [./temperature_gradient]
    type = ParsedFunction
    # T_top - scalar*y
    value = 10-(30)*y
  [../]

[]

[UserObjects]
active = 'initializeRBSystem performRBSystem'
#active = ''

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystemTransient
    parameters_filename = inputfiles/RB/ParallelLayers/TestModels/TestTransient.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    offline_stage = true
    execute_on = 'initial'
    #system = nl0
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineOnlineStageTransient

    exodus_file_name = TestTransient

    offline_stage = true
    online_stage = true
    offline_error_bound = false
    store_basis_functions = true

    mu_bar = 1
    online_mu = '1.05 2.5'

    execute_on = 'timestep_end'
    initial_rb_userobject = initializeRBSystem
    #system = nl0
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
parameter_names = 'mu_0 mu_1'

# Define the minimum and maximum value of the Theta object
mu_0 = '1.00 5.15'
mu_1 = '1.00 7.15'

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

# ======================= Transient RB system parameters =======================

# number of time steps
n_time_steps = 80

# size of time steps
delta_t = 0.04

# Generalized Euler method parameter in [0,1], euler_theta=1 implies backward Euler
euler_theta = 1

