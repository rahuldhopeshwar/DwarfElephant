[Mesh]
 file = meshs/CG/NoFaultModel/no_fault_model.msh
 partitioner = centroid
 centroid_partitioner_direction = y
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
  [./RBConduction0]
    type = DwarfElephantRBDiffusion
    block = 0
    ID_Aq = 0
    matrix_seperation_according_to_subdomains = false
  [../]
  [./RBConduction1]
    type = DwarfElephantRBDiffusion
    block = 1
    ID_Aq = 1
    matrix_seperation_according_to_subdomains = false
  [../]
  [./RBConduction2]
    type = DwarfElephantRBDiffusion
    block = 2
    ID_Aq = 2
    matrix_seperation_according_to_subdomains = false
  [../]
[]

[BCs]
[./RBtop]
  type = DwarfElephantRBDirichletBC
  boundary = 1
  value = 0.00
[../]

[./RBbottom]
  type = DwarfElephantRBNeumannBC
  boundary = 2
  value = 3.71
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
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg      101'
[]

[UserObjects]
[./initializeRBSystem]
  type = DwarfElephantInitializeRBSystemSteadyState
  execute_on = 'initial'
  N_max = 20
  n_training_samples = 100
  rel_training_tolerance = 1.e-7
  parameter_names = 'mu_0 mu_1 mu_2'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values = '0.50 0.50 0.50'
  parameter_max_values = '5.00 7.00 5.00'
[../]
[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  online_mu = '1.00 2.38 1.00'
  execute_on = 'timestep_end'
[../]
[]

[Outputs]
exodus = true
print_perf_log = true
execute_on = 'timestep_end'
  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]
