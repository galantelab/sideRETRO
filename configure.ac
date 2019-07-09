#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([sideRETRO], [0.4.0], [https://github.com/galantelab/sideRETRO/issues], [sider])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
AC_USE_SYSTEM_EXTENSIONS

LT_PREREQ([2.2])
LT_INIT([dlopen])

AM_INIT_AUTOMAKE

# Disable verbose gcc output
AM_SILENT_RULES([yes])

# Checks for programs.
AC_PROG_CC_C99
AC_PROG_INSTALL

# Checks for libraries.
PKG_CHECK_MODULES([CHECK], [check], ,
	[AC_MSG_WARN([check not found, It will not be possible to run 'make check'])])

AC_SEARCH_LIBS([dlopen], [dl], ,
	[AC_MSG_ERROR([libdl not found])])

AC_SEARCH_LIBS([gzopen], [z], ,
	[AC_MSG_ERROR([zlib not found. Please, manage to install it.
		For more details, see: https://github.com/madler/zlib])])

AC_SEARCH_LIBS([pthread_create], [pthread], ,
	[AC_MSG_ERROR([pthread (libpthread) not found.
		Please manage to install it])])

AC_SEARCH_LIBS([cos], [m], ,
	[AC_MSG_ERROR([math (libm) not found])])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_ASSERT
AC_CHECK_HEADERS([check.h zlib.h pthread.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Add third party library htslib and set the parameters
# passed to its configure script
AX_SUBDIRS_CONFIGURE([htslib],
	[[--disable-bz2], [--disable-lzma], [--disable-gcs], [--disable-libcurl], [--disable-s3]])

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT

echo \
"-------------------------------------------------

	${PACKAGE_NAME} Version ${PACKAGE_VERSION}

	Prefix: '${prefix}'.
	Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'
	Libraries: '${LIBS}'

	Now type 'make @<:@<target>@:>@'
	where the optional <target> is:
	all       - build all binaries
	install   - install everything

--------------------------------------------------"
