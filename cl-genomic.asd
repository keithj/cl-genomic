;;;
;;; Copyright (C) 2007-2009 Keith James. All rights reserved.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :cl-user)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (when (asdf:find-system :cl-system-utilities nil)
    (asdf:operate 'asdf:load-op :cl-system-utilities)))

(defpackage #:cl-genomic-system
  (:use :common-lisp :asdf :cl-system-utilities))

(in-package #:cl-genomic-system)

(defsystem cl-genomic
    :name "Common Lisp Genomics"
    :author "Keith James"
    :version "0.3.0"
    :licence "GPL v3"
    :depends-on (:trivial-gray-streams
                 :split-sequence
                 :cl-ppcre
                 :puri
                 :cl-gp-utilities
                 :cl-io-utilities)
    :in-order-to ((test-op (load-op :cl-genomic :cl-genomic-test)))
    :components ((:module :core
                          :pathname "src/"
                          :components
                          ((:file "package")
                           (:file "generics"
                            :depends-on ("package"))
                           (:file "conditions"
                            :depends-on ("package"))
                           (:file "bio-sequence-encoding"
                            :depends-on ("package"))
                           (:file "bio-alphabets"
                            :depends-on ("package"
                                         "generics"
                                         "bio-sequence-encoding"))
                           (:file "bio-sequence-classes"
                            :depends-on ("package"
                                         "bio-alphabets"
                                         "bio-sequence-encoding"))
                           (:file "bio-sequence-interval"
                            :depends-on ("package"
                                         "generics"
                                         "bio-sequence-classes"))
                           (:file "genetic-codes"
                            :depends-on ("package"
                                         "bio-alphabets"
                                         "bio-sequence-encoding"))
                           (:file "bio-sequence"
                            :depends-on ("package"
                                         "generics"
                                         "bio-alphabets"
                                         "bio-sequence-encoding"
                                         "bio-sequence-classes"
                                         "genetic-codes"))))
                 (:module :io
                          :pathname "src/io/"
                          :components
                          ((:file "bio-sequence-io-classes")
                           (:file "bio-sequence-io"
                            :depends-on ("bio-sequence-io-classes"))
                           (:file "fasta"
                            :depends-on ("bio-sequence-io"))
                           (:file "fastq"
                            :depends-on ("bio-sequence-io"))
                           (:file "format-conversion"
                            :depends-on ("bio-sequence-io"
                                         "fasta"
                                         "fastq"))
                           (:file "obo-io-classes")
                           (:file "obo-io"
                            :depends-on ("obo-io-classes"))
                           (:file "obo-powerloom"
                            :depends-on ("obo-io-classes"
                                         "obo-io"))
                           (:file "gff3"))
                          :depends-on (:core))
                 (:module :align
                          :pathname "src/alignment/"
                          :components
                          ((:file "bio-sequence-alignment")
                           (:file "matrices")                           
                           (:file "pairwise"
                            :depends-on ("matrices"
                                         "bio-sequence-alignment")))
                          :depends-on (:core))
                 (:lift-test-config :lift-tests
                                    :pathname "cl-genomic-test.config"
                                    :target-system :cl-genomic)
                 (:cldoc-config :cldoc-documentation
                                :pathname "doc/html"
                                :target-system :cl-genomic)))
