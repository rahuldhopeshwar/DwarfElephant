[Mesh]
 #file = RB_mesh_3layers.e
 type = GeneratedMesh
 dim = 3
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 zmin = 0
 zmax = 1
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
[./ RBtop]
  type = DwarfElephantRBDirichletBC
  boundary = 2
  value = 0.00
[../]

[./RBbottom]
  type = DwarfElephantRBNeumannBC
  boundary = 1
  value = 2.000000000000000e+01
[../]
[]

[Problem]
  type = DwarfElephantRBProblem
 # kernels = RBConduction
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
  parameter_names = 'mu_0'    #Please name them mu_0 , mu_1 , ..., mu_n for the reusability
  parameter_min_values = '1.0'
  parameter_max_values = '8.15'
[../]
[./ performRBSystem ]
  type = DwarfElephantOfflineOnlineStageSteadyState
  online_mu = '1.000000000000000e+00'
  execute_on = 'timestep_end'
[../]
[]

[Postprocessors]
  [./average]
    type = ElementAverageValue
    variable = temperature
    execute_on = 'custom'
  [../]
[]

[Outputs]
exodus = true
print_perf_log = true
  [./console]
    type = Console
    outlier_variable_norms = false
    output_postprocessors = false
  [../]
  [./DakotaOutput]
    type = DwarfElephantDakotaOutput
    file_path = '/home/dd823599/Dakota_first_example/'
    postprocessor = average
  [../]
[]
