ACLOCAL_AMFLAGS = -I m4

SUBDIRS = htslib bwa src tests

check-recursive:
	$(MAKE) -C tests check

README: README.md
	pandoc -f markdown -t plain --wrap=none $< -o $@

CLEANFILES = README
