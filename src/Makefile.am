bin_PROGRAMS = sider

sider_SOURCES = main.c

sider_LDFLAGS = -static

sider_LDADD = \
	./libsider.la \
	$(top_srcdir)/htslib/libhts.a \
	$(top_srcdir)/bwa/kopen.o \
	$(top_srcdir)/bwa/libbwa.a

noinst_LTLIBRARIES = libsider.la

libsider_la_SOURCES = \
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
		sam.h \
		bwa.c \
		bwa.h \
		align.c \
		align.h \
		sam_sort.c \
		sam_sort.h \
		abnormal.c \
		abnormal.h \
		bitree.c \
		bitree.h \
		ibitree.c \
		ibitree.h \
		gff.c \
		gff.h \
		str.c \
		str.h \
		sqlite3.c \
		sqlite3.h \
		db.c \
		db.h \
		chr.c \
		chr.h \
		exon.c \
		exon.h \
		thpool.c \
		thpool.h \
		process_sample.c \
		process_sample.h \
		io.c \
		io.h \
		logger.c \
		logger.h \
		dbscan.c \
		dbscan.h

libsider_la_CPPFLAGS = -I$(top_srcdir)/htslib
