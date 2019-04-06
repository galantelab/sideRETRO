TESTS = \
	check_sider \
	test_mem_leaks.sh

check_PROGRAMS = check_sider

EXTRA_DIST = test_mem_leaks.sh

check_sider_SOURCES = \
	check_sider.h \
	check_sider_main.c \
	check_sider_sider.c \
	check_sider_list.c \
	check_sider_hash.c \
	check_sider_array.c \
	check_sider_utils.c \
	check_sider_sam.c

check_sider_LDFLAGS = -static

check_sider_CFLAGS = $(CHECK_CFLAGS)

check_sider_LDADD = \
	$(top_builddir)/src/libsider.la \
	$(top_builddir)/htslib/libhts.a \
	$(CHECK_LIBS)
