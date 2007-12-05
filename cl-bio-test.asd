
(in-package #:cl-bio-system)

(defsystem cl-bio-test
    :depends-on (:cl-bio :fiveam)
    :components ((:module :cl-bio-test
                          :pathname "src/test/"
                          :components ((:file "package")
                                       (:file "bio-sequence-test"
                                              :depends-on ("package"))))))
