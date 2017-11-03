[Mesh]
 file = meshs/CG/FaultModelDislocated/fault_model_dislocated_refinement_0-12_0-06.e
[]

[Variables]
  [./temperature]
  [../]
[]

[GlobalParams]
  variable = temperature
[]

[Kernels]
  [./Conduction0]
    type = DwarfElephantFEThermalConduction
    thermal_conductivity = 1.00
    block = 0
  [../]
  [./Conduction1]
    type = DwarfElephantFEThermalConduction
    thermal_conductivity = 2.38
    block = 1
  [../]
  [./Conduction2]
    type = DwarfElephantFEThermalConduction
    thermal_conductivity = 1.00
    block = 2
  [../]
  [./ConductionFault2]
    type = DwarfElephantFEThermalConduction
    thermal_conductivity = 2.86
    block = 3
  [../]
  [./ConductionFault1]
    type = DwarfElephantFEThermalConduction
    thermal_conductivity = 2.86
    block = '4 5'
  [../]
[]

[BCs]
[./top]
  type = DirichletBC
  boundary = 2
  value = 0.00
[../]

[./bottom]
  type = NeumannBC
  boundary = 1
  value = 3.71
[../]
[]

[Executioner]
# Converges also without additional preconditioning settings
  type = Steady
  solve_type = 'Newton'
  l_tol = 1.0e-8
  nl_rel_tol = 1.0e-8
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg      101'
[]

[Outputs]
  exodus = true
  print_perf_log = true
  execute_on = 'timestep_end'
[]
