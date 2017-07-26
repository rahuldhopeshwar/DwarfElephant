[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
 # nz = 2
  xmin = 0.0
  xmax = 1
  ymin = 0.0
  ymax = 1
 # zmin = 0.0
 # zmax = 1
#  #elem_type=TET4
[]

[MeshModifiers]
#active = 'subdomains'
active = ' '
  [./subdomains]
    type = AssignElementSubdomainID
    subdomain_ids = '0 0 1 1'
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
active = 'RBConduction'
#active = 'Conduction'
#active = 'Conduction Euler'
  [./RBConduction]
    #type = DwarfElephantRBDiffusionLiftingFunction
    type = DwarfElephantRBDiffusionND
    variable = temperature
    initial_rb_userobject = initializeRBSystem
    u_ref = 10
    l_ref = 1000
    #lifting_function = temperature_gradient
    #simulation_type = transient
    #vector_seperation_according_to_subdomains = false
  [../]

  [./Conduction]
    #type = Conduction
    #type = Darcy
    type = DwarfElephantConductionLiftingFunction
    #type = DwarfElephantRBDiffusionLiftingFunction
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
#active = 'top'
#active = ' '
  [./RBtop]
    type = DwarfElephantRBFunctionDirichletBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 1 #4
    function = non_dimensionalize_temperature
    initial_rb_userobject = initializeRBSystem
    #simulation_type = transient
  [../]
  [./RBbottom]
    type = DwarfElephantRBNeumannBCND
    variable = temperature
    boundary = 3 #2
    value = -0.03
    u_ref = 10
    l_ref = 1000
    initial_rb_userobject = initializeRBSystem
    #simulation_type = transient
  [../]

  [./top]
    type = DirichletBC
    variable = temperature
    #boundary = 'lefttop righttop'
    boundary = 2
    value = 0
  [../]
  [./bottom]
    type = DirichletBC
    variable = temperature
    #boundary = 'leftbottom rightbottom'
    boundary = 0
    value = 0
  [../]
[]

[Materials]
active = ' '
#active = 'shale_top'
  [./shale_top]
    type = DwarfElephantShale
    block = 0
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
active = 'non_dimensionalize_temperature'
  [./temperature_gradient]
    type = ParsedFunction
    # T_top - scalar*y
    value = 10-(-30)*y
  [../]

  [./non_dimensionalize_temperature]
    type = ParsedFunction
    value = (10-10)/10
  [../]

[]

[UserObjects]
active = 'initializeRBSystem performRBSystem'
#active = ''

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystemSteadyState
    parameters_filename = inputfiles/RB/ParallelLayers/TestModels/TestNoDimensions.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    offline_stage = true
    execute_on = 'initial'
    #transient = true
    #system = nl0
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineOnlineStageSteadyState

    exodus_file_name = TestNoDimensions

    offline_stage = true
    online_stage = true
    offline_error_bound = false
    store_basis_functions = true

    mu_bar = 1
    online_mu = '1.05' # 2.5'

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
parameter_names = 'mu_0' # mu_1'

# Define the minimum and maximum value of the Theta object
mu_0 = '1.01 5.15'
#mu_1 = '1.01 7.15'

# Define the number of training sets for the Greedy-algorithm
n_training_samples = 10

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
n_time_steps = 10

# size of time steps
delta_t = 0.01

# Generalized Euler method parameter in [0,1], euler_theta=1 implies backward Euler
euler_theta = 1

