Introduction

cl-genomic is a portable Common Lisp library for processing genomic
data, including DNA, RNA and amino-acid sequences, and the biological
relationships between them.


Installation

cl-genomic uses ASDF for system definition. Copy or symlink
cl-genomic.asd (and optionally cl-genomic-test.asd) to your
asdf:*central-registry* and load cl-genomic with the asdf:operate
function:

 (asdf:operate 'asdf:load-op :cl-genomic)

or with the equivalent deoxybyte-systems:load-system function:

 (dxs:load-system :cl-genomic)


Tests

To run the unit and regression tests you need to have LIFT
installed. Run the tests with the asdf:operate function:

 (asdf:operate 'asdf:test-op :cl-genomic)

or with the equivalent deoxybyte-systems:test-system function:

 (dxs:test-system :cl-genomic)


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
                        git://git.b9.com/puri.git
CFFI                    http://common-lisp.net/project/cffi/
Ironclad                http://method-combination.net/lisp/ironclad/
                        git://github.com/froydnj/ironclad.git
cl-base64               http://www.cliki.net/cl-base64/
                        git://git.b9.com/cl-base64.git


Optional dependencies

LIFT                    http://common-lisp.net/project/lift/
                        git://github.com/gwkkwg/lift.git
CLDOC                   http://common-lisp.net/project/cldoc/
