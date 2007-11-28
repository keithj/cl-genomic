
(in-package :bigvector)

(eval-when (:compile-toplevel)
  (defconstant +bigvector-chunk-length+ (expt 2 19)
    "The length of the sequence of vectors which comprise a bigvector."))

(defclass bigvector ()
  ((chunks :initarg :chunks
           :reader chunks-of
           :documentation "A vector of chunks of equal size, except
for the last, which may be shorter than the others. Each chunk is
itself a specialised vector having the same element type."))
  (:documentation "A specialised vector of unlimited size."))

(defmacro create-chunks (bv-length element-type chunk-length)
  `(multiple-value-bind (full-chunks remainder)
    (floor ,bv-length ,chunk-length)
    (let ((chunks (loop for i from 0 below full-chunks
                        collect (make-array ,chunk-length
                                            :element-type ',element-type))))
      (unless (zerop remainder)
        (push (make-array remainder
                          :element-type ',element-type) chunks))
      (make-array (length chunks)
                  :element-type '(simple-array ,element-type ,chunk-length)
                  :initial-contents (nreverse chunks)))))

(defmethod length-of ((bv bigvector))
  (loop for chunk across (chunks-of bv)
        summing (length chunk) into bv-length
        finally (return bv-length)))

(defmethod bvref ((bv bigvector) (subscript integer))
  (address-cell (chunks-of bv) subscript))

(defun address-cell (data subscript)
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (declare (type simple-vector data)
           (type integer subscript))
  (multiple-value-bind (quotient remainder)
      (floor subscript +bigvector-chunk-length+)
    (let ((chunk (svref data quotient)))
      (declare (type (simple-array (unsigned-byte 2) (*)) chunk))
      (aref chunk remainder))))

(defmethod (setf bvref) (value (bv bigvector) (subscript integer))
  `(set-the-cell ,value (chunks-of ,bv) ,subscript))

(defun set-the-cell (value data subscript)
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (declare (type (unsigned-byte 2) value)
           (type simple-vector data)
           (type integer subscript))
  (multiple-value-bind (quotient remainder)
      (floor subscript +bigvector-chunk-length+)
    (let ((chunk (svref data quotient)))
      (declare (type (simple-array (unsigned-byte 2) (*)) chunk))
      (setf (aref chunk remainder) value))))

(defsetf address-cell (data subscript) (value)
  `(set-the-cell ,value ,data ,subscript))

(defmethod bigvector-subseq ((bv bigvector) start &optional
                             (end (length-of bv)))
  (cond ((< start 0)
         (error 'illegal-subscript :vector bv :subscript start))
        ((> end (length-of bv))
         (error 'subscript-out-of-bounds :vector bv :subscript end))
        ((< end start)
         (error 'incompatible-start-and-end :vector bv :start start :end end))
        (t
         (let ((sub (make-instance 'bigvector
                                   :length (- end start)
                                   :element-type (element-type-of bv))))
           (bigvector-copy bv sub :source-start start :source-end (1- end))))))

(defmethod bigvector-concatenate (&rest bigvectors)
  (let ((concat (make-instance 'bigvector
                               :length (reduce #'+  bigvectors :key #'length-of)
                               :element-type
                               (element-type-of (first bigvectors)))))
    (loop for bv in bigvectors
          for offset = 0 then (+ (length-of bv) offset)
          do (bigvector-copy bv concat :dest-start offset))
    concat))

(defmethod bigvector-copy ((source bigvector) (dest bigvector) &key
                           (source-start 0)
                           (source-end (1- (length-of source)))
                           (dest-start 0)
                           key)
  (loop for si from source-start to source-end
        for di = dest-start then (1+ di)
        do (setf (bvref dest di) (if key
                                     (funcall key (bvref source si))
                                   (bvref source si))))
  dest)

(defmethod bigvector-copy ((source bigvector) (dest vector) &key
                           (source-start 0)
                           (source-end (1- (length-of source)))
                           (dest-start 0)
                           key)
  (loop for si from source-start to source-end
        for di = dest-start then (1+ di)
        do (setf (aref dest di) (if key
                                    (funcall key (bvref source si))
                                  (bvref source si))))
  dest)

(defun array-copy (source target &key (source-start 0)
                   (source-end (1- (length source))) (target-start 0))
  (loop for si from source-start to source-end
        for ti = target-start then (1+ ti)
        do (setf (aref target ti) (aref source si)))
  target)
