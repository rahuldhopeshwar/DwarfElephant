[Mesh]
 file = RB_mesh_3layers.e
[]

[Variables]
  [./temperature]
  [../]
[]

[Kernels]
  [./RBConduction]
    type = DwarfElephantRBDiffusion
    variable = temperature
    initial_rb_userobject = initializeRBSystem
  [../]
[]

[BCs]
[./ RBtop]
  type = DwarfElephantRBDirichletBC
  variable = temperature
  boundary = 2
  value = 0.00
  initial_rb_userobject = initializeRBSystem
[../]

[./RBbottom]
  type = DwarfElephantRBNeumannBC
  variable = temperature
  boundary = 1
  value = 40
  initial_rb_userobject = initializeRBSystem
[../]
[]

[Problem]
  type = DwarfElephantRBProblem
[]

[Executioner]
  type = DwarfElephantRBExecutioner
  solve_type = 'Newton'
  l_tol = 1.0e-8
  nl_rel_tol = 1.0e-8
[]

[UserObjects]
[./initializeRBSystem]
  type = DwarfElephantInitializeRBSystemSteadyState
  execute_on = 'initial'
  N_max = 20
  n_training_samples = 100
  rel_training_tolerance = 1.e-5
  parameter_names = 'mu_0 mu_1 mu_2'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values = '1.0 1.0 1.0'
  parameter_max_values = '5.15 7.15 5.15'
[../]
[./ performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  store_basis_functions = true
  online_mu = '1.05 2.5 1.05'
  execute_on = 'timestep_end'
  initial_rb_userobject = initializeRBSystem
[../]
[]

[Outputs]
#print_perf_log = true
  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]
