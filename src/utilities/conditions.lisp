
(in-package :bigvector)

(define-condition illegal-subscript (error)
  ((vector :initarg :vector
           :reader vector-of)
   (subscript :initarg :subscript
              :reader subscript-of))
  (:report (lambda (condition stream)
              (format stream "~a is an illegal subscript for ~a."
                      (subscript-of condition)
                      (vector-of condition)))))

(define-condition incompatible-start-and-end (error)
  ((vector :initarg :vector
           :reader vector-of)
   (start :initarg :start
          :reader start-of)
   (end :initarg :end
        :reader end-of))
  (:report (lambda (condition stream)
             (format stream "End (~a) is smaller than start (~a) for ~a."
                     (start-of condition)
                     (end-of condition)
                     (vector-of condition)))))

(define-condition subscript-out-of-bounds (error)
  ((vector :initarg :vector
           :reader vector-of)
   (subscript :initarg :subscript
              :reader subscript-of))
  (:report (lambda (condition stream)
             (format stream "The subscript ~a exceeds the limit ~a for ~a."
                     (subscript-of condition)
                     (length (vector-of condition))
                     (vector-of condition)))))
