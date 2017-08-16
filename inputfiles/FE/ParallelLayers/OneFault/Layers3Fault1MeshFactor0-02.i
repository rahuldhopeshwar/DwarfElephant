[Mesh]
  file = meshs/parallel_3layers_1fault_meshfactor0-02.msh
[]

[MeshModifiers]
  [./AddFrac]
    type = MeshSideSet
    block_id = 3
    boundaries = 7
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
#active = 'RBConduction'
active = 'Conduction'
  [./RBConduction]
    type = DwarfElephantRBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
  [../]

  [./Conduction]
    type = DwarfElephantFEThermalConduction
    variable = temperature
  [../]
[]

[BCs]
#active = 'RBtop RBbottom'
active = 'top bottom'
  [./RBtop]
    type = DwarfElephantRBDirichletBC
    variable = temperature
    boundary = 4
    value = 0.00
    initial_rb_userobject = initializeRBSystem
    matrix_seperation_according_to_subdomains = true
  [../]

  [./RBbottom]
    type = DwarfElephantRBNeumannBC
    variable = temperature
    boundary = 2
    value = -40
    initial_rb_userobject = initializeRBSystem
    matrix_seperation_according_to_subdomains = true
  [../]

  [./top]
    type = DirichletBC
    variable = temperature
    boundary = 4
    value = 0
  [../]
  [./bottom]
    type = NeumannBC
    variable = temperature
    boundary = 2
    value = 40
  [../]

[Materials]
#active = ' '
active = 'shale_top sand_med shale_bot frac'
  [./shale_top]
    type = DwarfElephantShale
    block = 0
  [../]

   [./sand_med]
    type = DwarfElephantSandStone
    block = 1
  [../]

  [./shale_bot]
    type = DwarfElephantShale
    block = 2
  [../]

  [./frac]
    type = DwarfElephantFault
    block = 3
  [../]
[]

#[Problem]
#  type = DwarfElephantRBProblem
#[]

[Executioner]
  #type = DwarfElephantRBExecutioner
  type = Steady

  solve_type = 'Newton'
  l_tol = 1.e-8
  nl_rel_tol = 1.e-8
[]

[UserObjects]
#active = 'initializeRBSystem performRBSystem'
active = ''

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystemSteadyState
    parameters_filename = inputfiles/RB/ParallelLayers/SyntheticModels/OneFault/Layers3Fault1MeshFactor0-02.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    offline_stage = true
    execute_on = 'initial'
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineOnlineStageSteadyState

    exodus_file_name = Layers3Fault1MeshFactor0-02

    offline_stage = true
    online_stage = true
    offline_error_bound = false
    store_basis_functions = true

    mu_bar = 1
    online_mu = '1.05 2.5 1.05 3.0'

    execute_on = 'timestep_end'
    initial_rb_userobject = initializeRBSystem
  [../]
[]

[Outputs]
  print_perf_log = true
  exodus = true
  xda = false
  execute_on = 'timestep_end'

  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]


# ======================= Parameters for the RB method ========================

# Maximum number of basis functions that will be generated in the Offline-stage
Nmax = 50

# Name of the parameters
# Please name them mu_0, mu_1, ..., mu_n for the re-usability
parameter_names = 'mu_0 mu_1 mu_2 mu_3'

# Define the minimum and maximum value of the Theta object
mu_0 = '1.00 5.15'
mu_1 = '1.00 7.15'
mu_2 = '1.00 5.15'
mu_3 = '0.00 10.00'

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

