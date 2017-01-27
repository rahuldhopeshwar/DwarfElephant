[Mesh]
  file = syncline.msh
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]

  [./block1]
    type = Conduction
    variable = temperature
    block = '7 8 9'
  [../]

#  [./block2]
#    type = Conduction
#    variable = temperature
#    block = 8
#  [../]

#  [./block3]
#    type = Conduction
#    variable = temperature
#    block = 9
#  [../]
[]

[Materials]
  [./shale1]
    type = Shale
    block = 7
  [../]

  [./sandstone]
    type = SandStone
    block = 8
  [../]

  [./shale2]
    type = Shale
    block = 9
  [../]
[]

[BCs]
  [./bottom_shale1]
    type = DirichletBC
    variable = temperature
    boundary = 'bottom'
    value = 31
  [../]
  [./top_shale1]
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
  exodus = true
[]

