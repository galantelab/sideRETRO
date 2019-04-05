bin_PROGRAMS = sider
sider_SOURCES = main.c

sider_LDADD = ./libsider.la
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

libsider_la_LDFLAGS = -all-static
libsider_la_CPPFLAGS = -I$(top_srcdir)/htslib -L$(top_srcdir)/htslib
libsider_la_LIBADD = -lhts
