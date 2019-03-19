ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests

README: README.md
	pandoc -f markdown -t plain --wrap=none $< -o $@

CLEANFILES = README
