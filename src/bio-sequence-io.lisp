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

(in-package :bio-sequence)

(defmethod read-bio-sequence :before (stream format &key alphabet
                                             virtualp)
  (unless alphabet
    (error 'invalid-argument-error
           :params 'alphabet
           :args nil
           :text "an alphabet is required"))
  (unless (member virtualp '(t nil))
    (error 'invalid-argument-error
           :params 'virtualp
           :args virtualp
           :text "expected a value of T or NIL")))


(defun make-chunk-pname (file-pname chunk-number)
  "Returns a new pathname for a file chunk based on a file pathname
FILE-PNAME and an integer CHUNK-NUMBER."
  (make-pathname
   :directory (pathname-directory file-pname)
   :name (concatenate 'string (pathname-name file-pname) "."
                      (princ-to-string chunk-number))
   :type (pathname-type file-pname)))


(defmethod bio-sequence-io :before (format alphabet &optional handler
                                           &rest handler-initargs)
  (unless alphabet
    (error 'invalid-argument-error
           :params 'alphabet
           :args nil
           :text "an alphabet is required")))

(defgeneric begin-object (io-handler))
(defgeneric relation-property (io-handler relation value))
(defgeneric atomic-property (io-handler property value))
(defgeneric sequence-property (io-handler property index value))
(defgeneric end-object (io-handler))

(defclass io-handler ()
  ())

(defclass bio-sequence-handler (io-handler)
  ((atomic-props :initform nil
                 :accessor atomic-props-of)))

(defclass simple-sequence-handler (bio-sequence-handler)
  ((sequence-props :initform nil
                   :accessor sequence-props-of)))

(defclass quality-sequence-handler (simple-sequence-handler)
  ((metric :initarg :metric
           :reader metric-of)))

(defclass virtual-sequence-handler (bio-sequence-handler)
  ())

(defmethod begin-object ((handler bio-sequence-handler))
  (setf (atomic-props-of handler) nil
        (sequence-props-of handler) nil))

(defmethod begin-object ((handler virtual-sequence-handler)))

(defmethod end-object ((handler bio-sequence-handler)))

(defmethod atomic-property ((handler bio-sequence-handler)
                            (property symbol)
                            (value string))
  (setf (atomic-props-of handler)
        (acons property value (atomic-props-of handler))))

(defmethod atomic-property ((handler bio-sequence-handler)
                            (property symbol)
                            (value symbol))
  (setf (atomic-props-of handler)
        (acons property value (atomic-props-of handler))))

(defmethod sequence-property ((handler simple-sequence-handler)
                              (property symbol)
                              (index fixnum)
                              (value vector))
  (declare (ignore index))
  (let ((sequence-props (sequence-props-of handler)))
    (if (assoc property sequence-props)
        (vector-push-extend value (assocdr property sequence-props))
      (setf (sequence-props-of handler)
            (acons property
                   (make-array 1 :adjustable t :fill-pointer t
                               :initial-element value) sequence-props)))))

(defmethod sequence-property ((handler virtual-sequence-handler)
                              (property symbol)
                              (index fixnum)
                              (value vector))
  (when (eql :token-seq property)
    (let ((atomic-props (atomic-props-of handler)))
      (if (assoc :length atomic-props)
          (incf (assocdr :length atomic-props) (length value))
        (setf (atomic-props-of handler)
              (acons :length (length value) atomic-props))))))

(defmethod make-bio-sequence ((handler simple-sequence-handler))
  (let* ((atomic-props (atomic-props-of handler))
         (class (ecase (assocdr :alphabet atomic-props)
                  (:dna 'dna-sequence)
                  (:rna 'rna-sequence)))
         (sequence-props (sequence-props-of handler))
         (chunks (assocdr :token-seq sequence-props)))
    (when (zerop (length chunks))
      (error 'invalid-operation-error
             :text "attempt to make an empty concrete bio-sequence"))
    (let ((token-seq (etypecase (aref chunks 0)
                       (string (concat-strings chunks))
                       ((array (unsigned-byte 8))
                        (concat-into-sb-string chunks)))))
      (make-instance class
                     :identity (assocdr :identity atomic-props)
                     ;; :description (assocdr :description atomic-props)
                     :token-seq token-seq))))

(defmethod make-bio-sequence ((handler quality-sequence-handler))
  (let* ((atomic-props (atomic-props-of handler))
         (class (ecase (assocdr :alphabet atomic-props)
                  (:dna 'dna-quality-sequence)))
         (sequence-props (sequence-props-of handler))
         (sequence-chunks (assocdr :token-seq sequence-props))
         (quality-chunks (assocdr :quality sequence-props)))
    (when (zerop (length sequence-chunks))
      (error 'invalid-operation-error
             :text "no sequence data provided"))
    (when (zerop (length quality-chunks))
      (error 'invalid-operation-error
             :text "no quality data provided"))
    (let ((token-seq (etypecase (aref sequence-chunks 0)
                       (string (concat-strings sequence-chunks))
                       ((array (unsigned-byte 8))
                        (concat-into-sb-string sequence-chunks))))
          (quality (if (= 1 (length quality-chunks))
                       (aref quality-chunks 0)
                     (concat-quality-arrays quality-chunks))))
      (make-instance class
                     :identity (assocdr :identity atomic-props)
                     ;; :description (assocdr :description atomic-props)
                     :token-seq token-seq
                     :quality quality
                     :metric (metric-of handler)))))

(defmethod make-bio-sequence ((handler virtual-sequence-handler))
  (let* ((atomic-props (atomic-props-of handler))
         (class (ecase (assocdr :alphabet atomic-props)
                  (:dna 'dna-sequence)
                  (:rna 'rna-sequence))))
    (make-instance class
                   :identity (assocdr :identity atomic-props)
                   ;; :description (assocdr :description atomic-props)
                   :length (assocdr :length atomic-props))))

(defun concat-quality-arrays (quality-arrays)
  (let ((new-quality (make-array (reduce #'+ quality-arrays :key #'length)
                                 :element-type 'quality-score))
        (num-arrays (length quality-arrays)))
    (do ((i 0 (1+ i))
         (offset 0))
        ((= i num-arrays) new-quality)
      (let ((quality-array (aref quality-arrays i)))
        (unless (zerop (length quality-array))
          (copy-array quality-array 0 (1- (length quality-array))
                      new-quality offset)
          (incf offset (length quality-array)))))))
