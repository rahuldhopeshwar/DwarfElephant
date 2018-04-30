[Mesh]
 file = /Users/denise/projects/DwarfElephant/meshs/CG/NoFaultModel/no_fault_model_refinement_0-1.e
[]

[Variables]
  [./pressure]
  [../]
[]

[GlobalParams]
  initial_rb_userobject = initializeRBSystem
  variable = pressure
  simulation_type = transient
  min_x = 0.2
  max_x = 0.2
  min_y = 0.4
  max_y = 0.4
  min_z = 0.2
  max_z = 0.2
[]

[Kernels]
  [./RBConduction]
    type = DwarfElephantRBDiffusion
    compute_output = true
  [../]
  [./RBTime_0]
    type = DwarfElephantRBTimeDerivative
    compute_output = true
    ID_Mq = 1
  [../]
[]

[BCs]
[./RBleft]
  type = DwarfElephantRBDirichletBC
  boundary = 4
  value = 0.00
  ID_Mq = 1
  ID_Fq = 3
[../]

[./RBright]
  type = DwarfElephantRBPenaltyDirichletBC
  penalty = 100000
  boundary = 3
  value = 0.1
  split_boundary_according_to_subdomains = true
  subdomain_split = '0 1 2'
[../]
[]

[Problem]
  type = DwarfElephantRBProblem
[]

[Executioner]
  type = DwarfElephantRBExecutioner
  l_tol = 1.0e-10
  nl_rel_tol = 1.0e-10
  nl_rel_step_tol = 1.0e-10
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'
[]

[UserObjects]
[./initializeRBSystem]
  type = DwarfElephantInitializeRBSystemTransient
  execute_on = 'initial'
  N_max = 100
  n_training_samples = 100
  rel_training_tolerance = 1.e-7
  n_time_steps = 45
  delta_t = 1
  euler_theta = 1
  parameter_names = 'mu_0 mu_1 mu_2'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values = '0.5e-10 0.5e-08 0.5e-10'
  parameter_max_values = '1.5e-7 1.5e-5 1.5e-7'
[../]
[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageTransient
  online_stage = true
  online_mu = '1.e-3 1.e-1 1.e-3'
  execute_on = 'timestep_end'
  norm_online_values = true
  norm_id = 1
  compute_output = true
[../]
[]

[Outputs]
print_perf_log = true
execute_on = 'timestep_end'
exodus = false
  [./console]
    type = Console
    outlier_variable_norms = false
    execute_postprocessors_on = 'timestep_end'
  [../]
[]
