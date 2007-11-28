
(in-package :bigvector-test)

(define-test make-instance-bigvector-noargs
  "Test no-args construction."
  (let ((bv (make-instance 'bigvector)))
    (assert-eql 0 (length-of bv))
    (assert-eql t (element-type-of bv))))

(define-test make-instance-bigvector-args
  "Test with-args construction."
  (let* ((len 10)
         (type 'base-char)
         (bv (make-instance 'bigvector :length len :element-type type)))
    (assert-eql len (length-of bv))
    (assert-equal type (element-type-of bv))))

(define-test length-of-bigvector
  "Test length-of explicitly."
  (let* ((len 3500000)
         (type '(unsigned-byte 1))
         (bv (make-instance 'bigvector :length len :element-type type)))
    (assert-eql len (length-of bv))))

(define-test element-type-of-bigvector
  "Test element-type-of explicitly."
  (let* ((type '(unsigned-byte 1))
         (bv (make-instance 'bigvector :length 10 :element-type type)))
    (assert-equal type (element-type-of bv))))

(define-test bvref-bigvector
  "Test bvref, including setf."
  (let* ((len 3500000)
         (type '(unsigned-byte 1))
         (bv (make-instance 'bigvector :length len :element-type type)))
    (loop for i from 0 below (length-of bv)
          do (setf (bvref bv i) 1))
    (assert-true (count-elt 1 bv #'=))))

(define-test bigvector-subseq-bigvector
  "Test bigvector-subseq."
  (let* ((len 10)
         (type 'fixnum)
         (bv (make-instance 'bigvector :length len :element-type type)))
    (loop for i from 0 below (length-of bv)
          do (setf (bvref bv i) i))
    (assert-eql 'bigvector (type-of (bigvector-subseq bv 0)))

    (loop for i from 0 to (length-of bv)
          for sv = (bigvector-subseq bv 0 i)
          do (assert-eql i (length-of sv))
            (loop for j from 0 below (length-of sv)
                  do (assert-eql j (bvref sv j))))))

(define-test bigvector-concatenate-bigvector
  "Test bigvector-copy."
  (let* ((len 10)
         (type 'fixnum)
         (bv1 (make-instance 'bigvector :length len :element-type type))
         (bv2 (make-instance 'bigvector :length len :element-type type))
         (bv3 (make-instance 'bigvector :length len :element-type type)))
    (loop for i from 0 below len
          do (setf (bvref bv1 i) 1))
    (loop for i from 0 below len
          do (setf (bvref bv2 i) 2))
    (loop for i from 0 below len
          do (setf (bvref bv3 i) 3))

    (let ((concat (bigvector-concatenate bv1 bv2 bv3)))
      (loop for i from 0 below len
            do (assert-eql (bvref bv1 i) (bvref concat i)))
      (loop for i from 0 below len
            for j from (* 1 len) below (* 2 len)
            do (assert-eql (bvref bv2 i) (bvref concat j)))
      (loop for i from 0 below len
            for j from (* 2 len) below (* 3 len)
            do (assert-eql (bvref bv3 i) (bvref concat j))))))

(define-test bigvector-copy-bigvector
  "Test bigvector-copy."
  (let* ((len 10)
         (type 'fixnum)
         (bv1 (make-instance 'bigvector :length len :element-type type))
         (bv2 (make-instance 'bigvector :length len :element-type type))
         (v1 (make-array len :element-type type)))
    (loop for i from 0 to (1- (length-of bv1))
          do (setf (bvref bv1 i) i))

    ;; Copy to another bigvector
    (bigvector-copy bv1 bv2)
    (loop for i from 0 below (length-of bv1)
          do (assert-eql (bvref bv1 i) (bvref bv2 i)))

    ;; Copy to a regular array
    (bigvector-copy bv1 v1)
    (loop for i from 0 below (length-of bv1)
          do (assert-eql (bvref bv1 i) (aref v1 i)))))
    

(defun count-elt (x bigvector &optional (test #'eql))
  "Return the number of occurrences of X in BIGVECTOR using the
equality function TEST."
  (loop for i from 0 below (length-of bigvector)
        counting (funcall test x (bvref bigvector i)) into count
        finally return count))
