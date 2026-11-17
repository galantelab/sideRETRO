doc_sources = files(
  'authors.rst',
  'conf.py',
  'index.rst',
  'install.rst',
  'intro.rst',
  'method.rst',
  'output.rst',
  'result.rst',
  'retrocopy.rst',
  'usage.rst'
)

doc_images = files(
  'images/DBSCAN.svg',
  'images/LINE1.svg',
  'images/abnormal_alignment_chr.svg',
  'images/abnormal_alignment_clustering.svg',
  'images/abnormal_alignment_dist.svg',
  'images/abnormal_alignment_exon.svg',
  'images/abnormal_alignment_sr.svg',
  'images/alignment_of_CACNA1B.png',
  'images/abyzov_table2.png',
  'images/barbara.jpg',
  'images/breakpoints.svg',
  'images/intron_retention.svg',
  'images/polya_polarity.svg',
  'images/retrocopy_ends.svg',
  'images/genotype.svg',
  'images/indistinguishable_alignment.svg',
  'images/orientation_opposite_strand.svg',
  'images/orientation_same_strand.svg',
  'images/result_confusion.png',
  'images/retrocopy.svg'
)

sphinx_source_dir = join_paths(meson.global_source_root(), 'docs')
sphinx_builddir   = join_paths(meson.global_build_root(),  'docs')

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
           install_dir : join_paths(get_option('mandir'), 'man1'),
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
