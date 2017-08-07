[Mesh]
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
active = 'potential'
  [./potential]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
active = 'RBERT'
#active = 'ElectricalConduction0 ElectricalConduction1'
  [./RBERT]
    type = DwarfElephantRBDiffusion
    variable = potential
    initial_rb_userobject = initializeRBSystem
  [../]

  [./ElectricalConduction1]
    type = DwarfElephantFEElectricalConduction
    variable = potential
    resistivity = 200
    block = 1
  [../]

  [./ElectricalConduction0]
    type = DwarfElephantFEElectricalConduction
    variable = potential
    resistivity = 400
    block = 0
  [../]
[]

[BCs]
active = 'RBtop RBbottom'
#active = 'top bottom'
#active = ' '
  [./RBtop]
    type = DwarfElephantRBDirichletBC
    variable = potential
    #boundary = 'lefttop righttop'
    boundary = 3 #4
    value = 0.00
    initial_rb_userobject = initializeRBSystem
    ID_Aq = 1
  [../]
  [./RBbottom]
    type = DwarfElephantRBNeumannBC
    variable = potential
    boundary = 1 #2
    value = -1
    initial_rb_userobject = initializeRBSystem
    ID_Aq = 0
  [../]

  [./top]
    type = DirichletBC
    variable = potential
    boundary = 3
    value = 0
  [../]
  [./bottom]
    type = DirichletBC
    variable = potential
    boundary = 1
    value = 1
  [../]
[]

[Materials]
active = ' '
#active = 'shale_top sand_bottom'
  [./sand_bottom]
    type = DwarfElephantSandStone
    block = 0
  [../]

 [./shale_top]
   type = DwarfElephantShale
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

[UserObjects]
active = 'ERTPreCalculations initializeRBSystem performRBSystem'
#active = ''

 [./ERTPreCalculations]
    type = DwarfElephantERTPreCalculations
    execute_on = 'initial'
    position_A_electrode = '1. 2. 3. 4. 5. 6. 7.'
    position_B_electrode = '4. 5. 6. 7. 8. 9. 10.'
    position_M_electrode = '2. 3. 4. 5. 6. 7. 8.'
    position_N_electrode = '3. 4. 5. 6. 7. 8. 9.'
  [../]

  [./initializeRBSystem]
    type = DwarfElephantInitializeRBSystemSteadyState
    parameters_filename = inputfiles/RB/ParallelLayers/TestModels/TestERT.i
    skip_matrix_assembly_in_rb_system = true
    skip_vector_assembly_in_rb_system = true
    offline_stage = true
    execute_on = 'initial'
  [../]

  [./performRBSystem]
    type = DwarfElephantOfflineOnlineStageSteadyState

    exodus_file_name = TestERT

    offline_stage = true
    online_stage = true
    offline_error_bound = false
    store_basis_functions = true

    mu_bar = 1
    online_mu = '400 200'

    execute_on = 'timestep_end'
    initial_rb_userobject = initializeRBSystem
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
mu_0 = '380 420'
mu_1 = '180 220'

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

