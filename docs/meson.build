doc_sources = files(
  'authors.rst',
  'conf.py',
  'index.rst',
  'install.rst',
  'intro.rst',
  'method.rst',
  'result.rst',
  'retrocopy.rst',
  'usage.rst'
)

doc_images = files(
  'images/DBSCAN.png',
  'images/LINE1.png',
  'images/abnormal_alignment_chr.png',
  'images/abnormal_alignment_clustering.png',
  'images/abnormal_alignment_dist.png',
  'images/abnormal_alignment_exon.png',
  'images/abnormal_alignment_sr.png',
  'images/alignment_of_CACNA1B.png',
  'images/abyzov_table2.png',
  'images/barbara.jpg',
  'images/genotype.png',
  'images/indistinguishable_alignment.png',
  'images/logo_sideRETRO.png',
  'images/orientation_opposite_strand.png',
  'images/orientation_same_strand.png',
  'images/result_confusion.png',
  'images/retrocopy.png'
)

sphinx_source_dir = join_paths(meson.source_root(), 'docs')
sphinx_builddir   = join_paths(meson.build_root(),  'docs')

# Test if sphinx is installed
if sphinx_build.found()
  sphinx_target = get_option('sphinx_target')

  sphinx_cmd = [
    sphinx_build,
    '-b', sphinx_target,
    sphinx_source_dir,
    sphinx_builddir
  ]

  # If the user passed the option '-Dsphinx_target=man', then
  # we build a man page and set the installation to
  # /usr/local/man/man1
  if sphinx_target == 'man'
    message('Generate man page alongside the executable')

    custom_target('docs',
                 input : [doc_sources, doc_images],
                output : ['sider.1'],
               install : true,
      build_by_default : true,
           install_dir : join_paths('man', 'man1'),
               command : sphinx_cmd
    )
  else
    custom_target('docs',
                 input : [doc_sources, doc_images],
                output : ['.'],
               install : false,
      build_by_default : false,
               command : sphinx_cmd
    )
  endif
else
  # Make users aware that without sphinx the documentation
  # cannot be generated
  warning('Documentation is disabled without \'sphinx\'')
endif
