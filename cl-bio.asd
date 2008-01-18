
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
    :depends-on (:cl-gp-utilities :cl-io-utilities :split-sequence
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
                           (:file "bio-graph"
                                  :depends-on ("package"))
                           (:file "bio-ontology"
                                  :depends-on ("package"
                                               "bio-graph"))
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
