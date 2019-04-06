bin_PROGRAMS = sider

sider_SOURCES = main.c

sider_LDFLAGS = -static

sider_LDADD = \
	./libsider.la \
	$(top_srcdir)/htslib/libhts.a

noinst_LTLIBRARIES = libsider.la

libsider_la_SOURCES = \
		sider.c \
		sider.h \
		wrapper.c \
		wrapper.h \
		list.c \
		list.h \
		types.h \
		hash.c \
		hash.h \
		log.c \
		log.h \
		array.c \
		array.h \
		utils.c \
		utils.h \
		sam.c \
		sam.h

libsider_la_CPPFLAGS = -I$(top_srcdir)/htslib
