assembly_sources = files(
  'sucateiro.c'
)

assembler = 'sucateiro'

assembly_exec = executable(
  assembler,
  assembly_sources,
  include_directories : inc,
         dependencies : deps,
            link_with : sider_lib,
              install : false
)
