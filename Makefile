
# This Makefile is, for the most part, a wrapper around the Lisp build
# system ASDF. It provides targets for building non-Lisp components,
# such as texinfo documentation.


.PHONY:	all fasl doc coverage clean

default: all

all: fasl doc

fasl:
	sbcl --noinform --noprint \
	--eval "(progn (asdf:operate 'asdf:compile-op :cl-genomic) (quit))"

doc:
	sbcl --noinform --noprint --load make-doc.lisp

test:
	sbcl --noinform --noprint \
	--eval "(progn (asdf:operate 'asdf:test-op :cl-genomic) (quit))"

coverage: clean
	sbcl --noinform --noprint --load make-coverage-report.lisp

clean:
	find . -name \*.fasl -o -name \*.ufsl -o -name \*.lx64fsl -exec rm {} \;

