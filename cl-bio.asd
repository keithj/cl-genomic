
(defsystem cl-bio
    :name "Common Lisp Bioinformatics"
    :author "Keith James"
    :version "0.1"
    :licence "GPL"
    :depends-on (:cl-gp-utilities)
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
                                               "generics"))))))

