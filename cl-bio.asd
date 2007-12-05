
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
    :depends-on (:cl-gp-utilities :cl-io-utilities)
    :components ((:module :cl-bio
                          :pathname "src/"
                          :components
                          ((:file "package")
                           (:file "classes"
                                  :depends-on ("package"
                                               "bio-sequence-encoding"))
                           (:file "bio-sequence-encoding"
                                  :depends-on ("package"))
                           (:file "generics"
                                  :depends-on ("package"))
                           (:file "bio-sequence"
                                  :depends-on ("package"
                                               "classes"
                                               "bio-sequence-encoding"
                                               "generics"))
                           (:file "bio-sequence-io"
                                  :depends-on ("package"
                                               "classes"
                                               "bio-sequence-encoding"
                                               "generics"
                                               "bio-sequence"))))))

(in-package #:asdf)

(defmethod perform ((op test-op) (c (eql (find-system
                                          'cl-bio))))
  (operate 'load-op :cl-bio-test)
  (funcall (intern (string :run!) (string :fiveam))
           'cl-bio-system:testsuite))

(defmethod operation-done-p ((op test-op) (c (eql (find-system
                                                   'cl-bio))))
  nil)
