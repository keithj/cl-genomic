;;;
;;; Copyright (C) 2007-2008, Keith James. All rights reserved.
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

(defpackage #:cl-genomic-system
  (:use :common-lisp :asdf)
  (:export #:testsuite))


(in-package #:cl-genomic-system)

(defsystem cl-genomic
    :name "Common Lisp Bioinformatics"
    :author "Keith James"
    :version "0.1.1"
    :licence "GPL v3"
    :depends-on (:trivial-gray-streams
                 :split-sequence
                 :cl-ppcre
                 :puri
                 :cl-gp-utilities
                 :cl-io-utilities)
    :components ((:module :core
                          :pathname "src/"
                          :components
                          ((:file "package")
                           (:file "conditions"
                                  :depends-on ("package"))
                           (:file "bio-sequence-encoding"
                                  :depends-on ("package"))
                           (:file "bio-alphabets"
                                  :depends-on ("package"
                                               "bio-sequence-encoding"))
                           (:file "bio-sequence-classes"
                                  :depends-on ("package"
                                               "bio-alphabets"
                                               "bio-sequence-encoding"))
                           (:file "bio-sequence-interval"
                                  :depends-on ("package"
                                               "generics"
                                               "bio-sequence-classes"))
                           (:file "generics"
                                  :depends-on ("package"))
                           (:file "bio-sequence"
                                  :depends-on ("package"
                                               "generics"
                                               "bio-alphabets"
                                               "bio-sequence-encoding"
                                               "bio-sequence-classes"))))
                 (:module :io
                          :pathname "src/io/"
                          :components
                          ((:file "io-classes")
                           (:file "bio-sequence-io"
                                  :depends-on ("io-classes"))
                           (:file "fasta"
                                  :depends-on ("bio-sequence-io"))
                           (:file "fastq"
                                  :depends-on ("bio-sequence-io"))
                           (:file "format-conversion"
                                  :depends-on ("bio-sequence-io"
                                               "fasta"
                                               "fastq"))
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
                          :depends-on (:core))))


(in-package #:asdf)

(defmethod perform ((op test-op) (c (eql (find-system
                                          :cl-genomic))))
  (operate 'load-op :cl-genomic-test)

  (let ((*default-pathname-defaults* (component-pathname c)))
    (funcall (intern (string :run-tests) (string :lift))
             :config "cl-genomic-test.config")))

(defmethod operation-done-p ((op test-op) (c (eql (find-system
                                                   :cl-genomic))))
  nil)

(defmethod perform ((op cldoc-op) (c (eql (find-system
                                           :cl-genomic))))
  (unless (find-package :bio-sequence)
    (operate 'load-op :cl-genomic))

  (let ((*default-pathname-defaults* (component-pathname c))
        (fn-sym (intern (string :extract-documentation) (string :cldoc)))
        (op-sym (intern (string :html) (string :cldoc))))
    (funcall fn-sym op-sym "./doc/html" c)))
