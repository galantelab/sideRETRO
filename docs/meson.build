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
  'images/barbara.jpg',
  'images/genotype.png',
  'images/indistinguishable_alignment.png',
  'images/logo_sideRETRO.png',
  'images/orientation_opposite_strand.png',
  'images/orientation_same_strand.png',
  'images/result_confusion.png',
  'images/result_heatmap.png',
  'images/retrocopy.png'
)

sphinx_source_dir = join_paths(meson.source_root(), 'docs')
sphinx_builddir   = join_paths(meson.build_root(), 'meson-docs')

# Test if sphinx is installed
if sphinx_build.found()
  custom_target('docs',
      input : [doc_sources, doc_images],
     output : ['.'],
    command : [sphinx_build, '-M', get_option('sphinx_target'),
      sphinx_source_dir, sphinx_builddir]
  )
else
  # Make users aware that without sphinx the documentation
  # cannot be generated
  warning('Documentation is disabled without \'sphinx\'')
endif
