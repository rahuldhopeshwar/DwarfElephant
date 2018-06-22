[Mesh]
 file = RB_mesh_3layers.e
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

[BCs]
[./RBtop]
  type = DwarfElephantRBDirichletBC
  boundary = 2
  value = 0.00
[../]

[./RBbottom]
  type = DwarfElephantRBNeumannBC
  boundary = 1
  value = 40
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
  normalize_rb_bound_in_greedy = true
[../]
[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  online_mu = '1.05 2.5 1.05'
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
