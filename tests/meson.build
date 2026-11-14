test_sources = files(
  'check_sider.h',
  'check_sider_abnormal.c',
  'check_sider_array.c',
  'check_sider_bed.c',
  'check_sider_bitree.c',
  'check_sider_blacklist.c',
  'check_sider_chr.c',
  'check_sider_cigar.c',
  'check_sider_cluster.c',
  'check_sider_correlation.c',
  'check_sider_db.c',
  'check_sider_db_merge.c',
  'check_sider_dbscan.c',
  'check_sider_dedup.c',
  'check_sider_exon.c',
  'check_sider_fault.c',
  'check_sider_fasta.c',
  'check_sider_genotype.c',
  'check_sider_gff.c',
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
  'check_sider_summary.c',
  'check_sider_utils.c',
  'check_sider_vcf.c',
  'check_sider_wrapper.c'
)

# Package executable tester
check = '_'.join(['check', package])

# Fault injection: Wrap the allocators and the pipes at the
# check executable only, if the linker supports it (GNU ld, lld).
# The wrappers are at check_sider_fault.c
wrap_args = []
wrap_c_args = []

foreach f : ['malloc', 'calloc', 'realloc', 'strdup', 'popen', 'pclose',
  'regcomp', 'pthread_mutex_init', 'gzclose', 'bam_init1', 'bam_dup1',
  'bam_copy1', 'hts_close', 'sam_itr_queryi', 'hts_set_threads',
  'sam_hdr_write', 'sam_read1', 'hdopen', 'hts_hopen']
  wrap_args += '-Wl,--wrap=' + f
endforeach

if cc.has_multi_link_arguments(wrap_args)
  wrap_c_args += '-DHAVE_LD_WRAP'
else
  wrap_args = []
  warning('Out of memory tests are disabled: The linker has no \'--wrap\'')
endif

# Test only if libcheck is found
if check_dep.found()
  # The executable tester for the units
  e = executable(
    check,
    test_sources,
    include_directories : inc,
           dependencies : [deps, check_dep],
           link_with    : sider_lib,
                 c_args : wrap_c_args,
              link_args : wrap_args,
                install : false
  )

  # Enable testing
  test('all', e, suite : 'unit')
else
  # Make users aware that without libcheck the test
  # does not run
  warning('Testing is disabled without \'libcheck\'')
endif

# Simulated retrocopies through the whole pipeline:
# process-sample, merge-call and make-vcf
python3 = find_program('python3', required : false)

if python3.found()
  test(
    'simulation',
    python3,
       args : [
         files('sim/run_sim_test.py'),
         sider_exec.full_path(),
         files('sim/retrocopy_sim.py')
       ],
    depends : sider_exec,
      suite : 'simulation',
    timeout : 600
  )

  test(
    'cli',
    python3,
       args : [
         files('cli/run_cli_test.py'),
         sider_exec.full_path(),
         files('sim/retrocopy_sim.py')
       ],
    depends : sider_exec,
      suite : 'cli',
    timeout : 600
  )
else
  warning('Simulation and CLI tests are disabled without \'python3\'')
endif
