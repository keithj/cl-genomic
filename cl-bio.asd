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

(defpackage #:cl-bio-system
  (:use :common-lisp :asdf)
  (:export #:testsuite))


(in-package #:cl-bio-system)

(defsystem cl-bio
    :name "Common Lisp Bioinformatics"
    :author "Keith James"
    :version "0.1.1"
    :licence "GPL v3"
    :depends-on (:cl-gp-utilities :cl-io-utilities
                                  :trivial-gray-streams :split-sequence
                                  :cl-ppcre :puri)
    :components ((:module :core
                          :pathname "src/"
                          :components
                          ((:file "package")
                           (:file "bio-sequence-encoding"
                                  :depends-on ("package"))
                           (:file "classes"
                                  :depends-on ("package"
                                               "bio-sequence-encoding"))
                           (:file "generics"
                                  :depends-on ("package"))
                           (:file "frames"
                                  :depends-on ("package"))
                           (:file "bio-sequence"
                                  :depends-on ("package"
                                               "generics"
                                               "bio-sequence-encoding"
                                               "classes"))
                           (:file "bio-sequence-io"
                                  :depends-on ("package"))))
                 (:module :io
                          :pathname "src/io/"
                          :components
                          ((:file "fasta")
                           (:file "fastq")
                           (:file "gff3"))
                          :depends-on (:core))))


(in-package #:asdf)

(defmethod perform ((op test-op) (c (eql (find-system
                                          'cl-bio))))
  (operate 'load-op :cl-bio-test)
  (let ((*default-pathname-defaults* (component-pathname c)))
    (funcall (intern (string :run!) (string :fiveam))
             'cl-bio-system:testsuite)))

(defmethod operation-done-p ((op test-op) (c (eql (find-system
                                                   'cl-bio))))
  nil)
