[Mesh]
 file = libMesh_EIM_example_mesh.msh
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
[]

[KernelsEIMFAction]
  #variable = temperature
  initialize_rb_system = initializeRBSystem
[]

[BCs] # To be changed
[./top]
  type = DwarfElephantRBDirichletBC
  boundary = 1
  value = 0.00
[../]

[./bottom]
  type = DwarfElephantRBNeumannBC
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
[]

[UserObjects] # for user objects with the same execute_on option, the order of execution will be decided by the alphabetical order of the userobject sub-block names. i. e. in this case, the object initializeRBSystem will be executed before the object jEIMInnerProductMatrixComputation

[./initializeRBSystem]
  type = DwarfElephantInitializeRBSystemSteadyState
  execute_on = 'initial'
  N_max_EIM = 20
  n_training_samples_EIM = 25
  rel_training_tolerance_EIM = 0.001
  parameter_names_EIM = 'center_x center_y'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values_EIM = '-1.0 -1.0'
  parameter_max_values_EIM = '1.0 1.0'
  best_fit_type_EIM = projection
[../]

[./jEIMInnerProductMatrixComputation]
  type = DwarfElephantComputeEIMInnerProductMatrixSteadyState
  execute_on = 'initial'
  n_training_samples_RB = 100
  N_max_RB = 15
  parameter_names_RB = 'center_x center_y'
  parameter_min_values_RB = '-1.0 -1.0'
  parameter_max_values_RB = '1.0 1.0'
  deterministic_training_RB = true 
  initialize_rb_userobject = initializeRBSystem
[../]

[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  online_mu = '0.3 -0.4'
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
