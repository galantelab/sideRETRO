test_sources = files(
  'check_sider.h',
  'check_sider_abnormal.c',
  'check_sider_array.c',
  'check_sider_bed.c',
  'check_sider_bitree.c',
  'check_sider_blacklist.c',
  'check_sider_chr.c',
  'check_sider_cluster.c',
  'check_sider_correlation.c',
  'check_sider_db.c',
  'check_sider_db_merge.c',
  'check_sider_dbscan.c',
  'check_sider_dedup.c',
  'check_sider_exon.c',
  'check_sider_fasta.c',
  'check_sider_genotype.c',
  'check_sider_gff.c',
  'check_sider_graph.c',
  'check_sider_gz.c',
  'check_sider_hash.c',
  'check_sider_ibitree.c',
  'check_sider_io.c',
  'check_sider_list.c',
  'check_sider_main.c',
  'check_sider_process_sample.c',
  'check_sider_retrocopy.c',
  'check_sider_sam.c',
  'check_sider_set.c',
  'check_sider_set.c',
  'check_sider_str.c',
  'check_sider_utils.c',
  'check_sider_vcf.c',
  'check_sider_wrapper.c'
)

# Package executable tester
check = '_'.join(['check', package])

# Test only if libcheck is found
if check_dep.found()
  # The executable tester for the units
  e = executable(
    check,
    test_sources,
    include_directories : inc,
           dependencies : [deps, check_dep],
           link_with    : sider_lib,
                install : false
  )

  # Enable testing
  test('all', e)
else
  # Make users aware that without libcheck the test
  # does not run
  warning('Testing is disabled without \'libcheck\'')
endif
