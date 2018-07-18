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
  N_max = 20 # in libmesh input file
  n_training_samples = 100 # in libmesh input file
  rel_training_tolerance = 1.e-5 # not used in libmesh input file
  parameter_names = 'mu_0 mu_1 mu_2'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability (not used in libmesh input file)
  parameter_min_values = '1.0 1.0 1.0' # used in libmesh input file
  parameter_max_values = '5.15 7.15 5.15' # used in libmesh input file
[../]
[./example_gen_usr_obj]
  type = ExampleGeneralUserObject
  variable = temperature
  execute_on = 'initial'
[../]
[./performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  online_mu = '1.05 2.5 1.05' # used in libmesh input file
  # parameter online_N is used in libmesh input file, but not here
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
print_perf_log = false
  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]
