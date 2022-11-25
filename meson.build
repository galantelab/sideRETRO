project(
  'sideRETRO', 'c',
          version : '1.1.0',
          license : 'GPL3',
    meson_version : '>= 0.48.0',
  default_options : ['c_std=c99']
)

# Package executable
package = 'sider'

# Package name as in project
package_name = meson.project_name()

# Version template
version_t = '@version@'

# Package string as in autoheader
package_string = ' '.join([meson.project_name(), version_t])

# Package bugreport
package_bugreport = 'https://github.com/galantelab/sideRETRO/issues'

# Create a configuration structure
config = configuration_data()

# Defines
set_defines = [
  ['PACKAGE',           package],
  ['PACKAGE_NAME',      package_name],
  ['PACKAGE_BUGREPORT', package_bugreport],
  ['PACKAGE_STRING',    package_string],
  ['PACKAGE_VERSION',   version_t],
  ['VERSION',           version_t]
]

# Set config values
foreach define: set_defines
  config.set_quoted(define[0], define[1])
endforeach

# Write macros to 'config.h.in' file
config_in = configure_file(
         output : 'config.h.in',
  configuration : config
)

# The version is dynamically generated according to
# the 'git tag and commit'
# Read config.h.in and replace the version template
# generating in turn config.h
config_out = vcs_tag(
         command : ['git', 'describe', '--tags', '--dirty=+'],
           input : config_in,
          output : 'config.h',
  replace_string : version_t
)

# Declare 'config.h' as dependency in order to obligate
# the build system to generate it before compiling
# Fix issue #1
config_dep = declare_dependency(sources : config_out)

# Compiler flags
common_flags = ['-DHAVE_CONFIG_H', '-D_GNU_SOURCE']

# Set project flags and language
add_project_arguments(common_flags, language: 'c')

# Let the compiler find required libraries
cc = meson.get_compiler('c')

# Math library portably
m_dep = cc.find_library('m', required : false)

# Libdl is required by SQLite3
dl_dep  = cc.find_library('dl', required : false)

# Enable threads in a portable way. Let meson do the hard work
thread_dep = dependency('threads')

# The following dependencies, if not found in the host, will be
# downloaded and statically linked

# Zlib
zlib_dep = dependency('zlib', fallback : ['zlib', 'zlib_dep'])

# SQLite3
sqlite_dep = dependency('sqlite3', version : '>= 3.28.0', fallback : ['sqlite', 'sqlite_dep'])

# HTSLib
htslib_dep = dependency('htslib', version : '>= 1.9', fallback : ['htslib', 'htslib_dep'])

# Libcheck framework is used for testing, so it is not
# required for building and installing
check_dep = dependency('check', required : false)

# Sphinx documentation generator
sphinx_build = find_program('sphinx-build', required : false)

# All required deps
deps = [
  config_dep,
  thread_dep,
  sqlite_dep,
  htslib_dep,
  zlib_dep,
  m_dep,
  dl_dep
]

# Set the path to find 'config.h'
inc = include_directories('.')

# Build static library and executable
subdir('src')

# Testing
subdir('tests')

# Documentation
subdir('docs')
