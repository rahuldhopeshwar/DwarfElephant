[Mesh]
 file = libMesh_EIM_example_mesh2.msh
[]

[Variables]
  [./temperature]
  [../]
[]

[GlobalParams]
  initial_rb_userobject = initializeRBSystem
  variable = temperature
[]

[Kernels]
  [./RBConduction]
    type = DwarfElephantRBDiffusion
  [../]
  
  [./EIMF]
    type = DwarfElephantEIMFKernel
  [../]
[]

[BCs] # To be changed
[./top]
  type = DwarfElephantRBDirichletBC
  boundary = 1
  value = 0.00
[../]

[./bottom]
  type = DwarfElephantRBDirichletBC
  boundary = 3
  value = 0
[../]

[./left]
  type = DwarfElephantRBDirichletBC
  boundary = 2
  value = 0.00
[../]

[./right]
  type = DwarfElephantRBDirichletBC
  boundary = 4
  value = 0.00
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
  offline_stage = false
[]

[UserObjects]

[./initializeRBSystem]
  type = DwarfElephantInitializeRBSystemSteadyState
  use_EIM = true
  execute_on = 'initial'
  N_max_EIM = 20
  n_training_samples_EIM = 25
  rel_training_tolerance_EIM = 0.001
  parameter_names_EIM = 'mu_0 mu_1'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values_EIM = '-1.0 -1.0'
  parameter_max_values_EIM = '1.0 1.0'
  deterministic_training_EIM = true
  best_fit_type_EIM = projection
  n_training_samples_RB = 100
  N_max_RB = 9
  rel_training_tolerance_RB = 0.001
  parameter_names_RB = 'mu_0 mu_1'
  parameter_min_values_RB = '-1.0 -1.0'
  parameter_max_values_RB = '1.0 1.0'
  deterministic_training_RB = true
  offline_stage = false
[../]

[./jEIMInnerProductMatrixComputation]
  type = DwarfElephantComputeEIMInnerProductMatrixSteadyState
  execute_on = "EIM"
  initialize_rb_userobject = initializeRBSystem
[../]

[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  online_stage = true
  online_N = 9
  online_mu = '0.7 -0.3'
  execute_on = 'timestep_end'
[../]
[]

#[Postprocessors]
#  [./average]
#    type = ElementAverageValue
#    variable = temperature
#    execute_on = 'custom'
#  [../]
#[]

[Outputs]
exodus = true
# csv = true   # only required for the PostProcessors
print_perf_log = true
  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]
