TESTS = \
	check_sider \
	test_mem_leaks.sh

check_PROGRAMS = check_sider

EXTRA_DIST = test_mem_leaks.sh

check_sider_SOURCES = \
	check_sider.h \
	check_sider_main.c \
	check_sider_list.c \
	check_sider_hash.c \
	check_sider_array.c \
	check_sider_utils.c \
	check_sider_sam.c \
	check_sider_bwa.c \
	check_sider_bitree.c \
	check_sider_ibitree.c \
	check_sider_str.c \
	check_sider_db.c \
	check_sider_chr.c \
	check_sider_exon.c \
	check_sider_abnormal.c \
	check_sider_gff.c \
	check_sider_io.c \
	check_sider_process_sample.c

check_sider_LDFLAGS = -static

check_sider_CFLAGS = $(CHECK_CFLAGS)

check_sider_CPPFLAGS = -I$(top_srcdir)/htslib

check_sider_LDADD = \
	$(top_builddir)/src/libsider.la \
	$(top_builddir)/htslib/libhts.a \
	$(top_builddir)/bwa/kopen.o \
	$(top_builddir)/bwa/libbwa.a \
	$(CHECK_LIBS)
