
(defsystem cl-bio-test
    :depends-on (:cl-bio :lisp-unit)
    :components ((:module :bigvector-test
                          :pathname "src/test/"
                          :components ((:file "package")
                                       (:file "bio-sequence-test"
                                              :depends-on ("package")))))
    :in-order-to ((test-op (load-op cl-bio-test))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (unless (find-package 'lisp-unit)
    (asdf:oos 'asdf:load-op 'lisp-unit))
  (unless (find-package 'cl-bio)
    (asdf:oos 'asdf:load-op 'cl-bio)))

(defmethod perform ((operation test-op)
                    (component (eql (find-system 'cl-bio-test))))
  (lisp-unit:run-all-tests :cl-bio-test))

(defmethod operation-done-p ((operation test-op)
                             (component (eql (find-system 'cl-bio-test))))
  (values nil))
