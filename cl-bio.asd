
(in-package :cl-user)

(defpackage #:cl-bio-system
  (:use :common-lisp :asdf)
  (:export #:testsuite))


(in-package #:cl-bio-system)

(defsystem cl-bio
    :name "Common Lisp Bioinformatics"
    :author "Keith James"
    :version "0.1"
    :licence "GPL"
    :depends-on (:cl-gp-utilities :cl-io-utilities :split-sequence)
    :components ((:module :cl-bio
                          :serial t
                          :pathname "src/"
                          :components
                          ((:file "package")
                           (:file "bio-sequence-encoding")
                           (:file "classes")
                           (:file "generics")
                           (:file "bio-graph")
                           (:file "bio-sequence")
                           (:file "bio-sequence-io")
                           (:file "fasta")
                           (:file "fastq")))))


(in-package #:asdf)

(defmethod perform ((op test-op) (c (eql (find-system
                                          'cl-bio))))
  (operate 'load-op :cl-bio-test)
  (funcall (intern (string :run!) (string :fiveam))
           'cl-bio-system:testsuite))

(defmethod operation-done-p ((op test-op) (c (eql (find-system
                                                   'cl-bio))))
  nil)
