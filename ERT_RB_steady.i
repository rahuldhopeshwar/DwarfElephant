[Mesh]
 file = RB_mesh_3layers.e
[]

[Variables]
  [./potential]
  [../]
[]

[GlobalParams]
  initial_rb_userobject = initializeRBSystem
  variable = potential
[]

[Kernels]
  [./block_012]
    type = DwarfElephantRBDiffusion
  [../]
[]

[DiracKernels]
  [./source]
    type = DwarfElephantRBConstantPointSource
    matrix_seperation_according_to_subdomains = false
    ID_Aq = 3
    value = 1.e-3
    point = '0.3 0 0'
  [../]
  [./sink]
    type = DwarfElephantRBConstantPointSource
    matrix_seperation_according_to_subdomains = false
    ID_Aq = 3
    value = -1.e-3
    point = '0.6 0 0'
  [../]
[]

[BCs]
[./top]
  type = DwarfElephantRBPenaltyDirichletBC
  penalty = 100000
  boundary = 2
  value = 0.0
  ID_Fq = 1
[../]
[./bottom]
  type = DwarfElephantRBNeumannBC
  boundary = 1
  value = 0.1
[../]
[]

[Problem]
  type = DwarfElephantRBProblem
[]

[Executioner]
  type = DwarfElephantRBExecutioner
  solve_type = 'Newton'
  l_tol = 1.0e-10
  nl_rel_tol = 1.0e-10
[]

[UserObjects]
[./initializeRBSystem]
  type = DwarfElephantInitializeRBSystemSteadyState
  execute_on = 'initial'
  offline_stage = true
  N_max = 20
  n_training_samples = 100
  rel_training_tolerance = 1.e-7
  parameter_names = 'mu_0 mu_1 mu_2'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values = '0.50 0.10 0.50'
  parameter_max_values = '2.00 3.50 2.00'
[../]
[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  offline_stage = true
  online_stage = true
  online_mu = '1.00 2.00 1.00'
  execute_on = 'timestep_end'
  online_N = 12
  n_outputs = 1
[../]
[]

[Outputs]
execute_on = 'timestep_end'
exodus = true
  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]
