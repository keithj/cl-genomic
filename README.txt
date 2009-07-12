Introduction

cl-genomic is a portable Common Lisp library for processing genomic
data, including DNA, RNA and amino-acid sequences, and the biological
relationships between them.


Documentation

See the Lisp docstrings, particularly the package docstrings for an
overview. HTML documentation may be generated with the command:

 (dxs:document-system :cl-genomic)

at the REPL, provided that CLDOC is installed.

The manual provides an overview of the design philosophy and examples
of use.


Dependencies

deoxybyte-systems       git://github.com/keithj/deoxybyte-systems.git
deoxybyte-utilities     git://github.com/keithj/deoxybyte-utilities.git
deoxybyte-io            git://github.com/keithj/deoxybyte-io.git
deoxybyte-unix          git://github.com/keithj/deoxybyte-unix.git

cl-ppcre                http://weitz.de/cl-ppcre/
PURI                    http://puri.b9.com/
CFFI                    http://common-lisp.net/project/cffi/
