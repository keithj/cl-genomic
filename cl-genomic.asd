;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
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

(asdf:load-system :deoxybyte-systems)

(in-package :uk.co.deoxybyte-systems)

(defsystem cl-genomic
    :name "Common Lisp Genomics"
    :version "0.6.1"
    :author "Keith James"
    :licence "GPL v3"
    :depends-on ((:version :cffi "0.10.3")
                 (:version :cl-ppcre "2.0.0")
                 (:version :cl-base64 "3.1") ; Actually 3.3.1, but
                                             ; cl-base64's ASDF file
                                             ; is out of sync with its
                                             ; professed version
                 (:version :ironclad "0.27")
                 (:version :deoxybyte-io "0.6.2")
                 (:version :deoxybyte-unix "0.6.5")
                 (:version :deoxybyte-utilities "0.7.0")
                 :puri)
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
                           (:file "checksum"
                            :depends-on ("package"
                                         "bio-alphabets"
                                         "bio-sequence-encoding"
                                         "bio-sequence-classes"))
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
                                         "genetic-codes"))
                           (:file "fjoin"
                            :depends-on ("package"
                                         "bio-sequence-interval"))))
                 (:module :io
                          :pathname "src/io/"
                          :components
                          ((:file "bio-sequence-io-classes")
                           (:file "bio-sequence-io"
                            :depends-on ("bio-sequence-io-classes"))
                           (:file "pure"
                            :depends-on ("bio-sequence-io"))
                           (:file "raw"
                            :depends-on ("bio-sequence-io"))
                           (:file "fasta"
                            :depends-on ("bio-sequence-io"))
                           (:file "fastq"
                            :depends-on ("bio-sequence-io"))
                           (:file "format-conversion"
                            :depends-on ("bio-sequence-io"
                                         "pure"
                                         "raw"
                                         "fasta"
                                         "fastq"))
                           (:file "obo-io-classes")
                           (:file "obo"
                            :depends-on ("obo-io-classes"))
                           (:file "obo-powerloom"
                            :depends-on ("obo-io-classes"
                                         "obo"))
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
                                    :pathname "cl-genomic-test"
                                    :target-system :cl-genomic)
                 (:cldoc-config :cldoc-documentation
                                :pathname "doc/html/"
                                :target-system :cl-genomic)))
