doc_sources = files(
  'authors.rst',
  'index.rst',
  'install.rst',
  'intro.rst',
  'conf.py',
  'method.rst',
  'usage.rst'
)

sphinx_source_dir = join_paths(meson.source_root(), 'docs')
sphinx_builddir   = join_paths(meson.build_root(), 'meson-docs')

# Test if sphinx is installed
if sphinx_build.found()
  custom_target('docs',
      input : doc_sources,
     output : ['.'],
    command : [sphinx_build, '-M', get_option('sphinx_target'),
      sphinx_source_dir, sphinx_builddir]
  )
else
  # Make users aware that without sphinx the documentation
  # cannot be generated
  warning('Documentation is disabled without \'sphinx\'')
endif
