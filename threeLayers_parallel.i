[Mesh]
  file = parallel.msh
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]

  [./blocks]
    type = Conduction
    variable = temperature
  [../]
[]

[Materials]
  [./shale_top]
    type = Shale
    block = 7
  [../]

  [./sandstone]
    type = SandStone
    block = 8
  [../]

  [./shale_bottom]
    type = Shale
    block = 9
  [../]
[]

[BCs]
  [./bottom]
    type = DirichletBC
    variable = temperature
    boundary = 'bottom'
    value = 31
  [../]
  [./top]
    type = DirichletBC
    variable = temperature
    boundary = 'top'
    value = 10
  [../]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]

[Outputs]
  execute_on = timestep_end
  file_base = threeLayers_parallel_shale_1.05_sand_2.5
  exodus = true
  xda = true
[]

