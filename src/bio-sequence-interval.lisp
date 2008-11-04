;;;
;;; Copyright (C) 2008 Keith. All rights reserved.
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

;;; For practicality we allow the reference to be null so that we can
;;; have unattached intervals. An alternative is to use proxies, which
;;; may not be practical when there are millions of intervals and
;;; millions of references.

(defclass interval ()
  ((reference :initform nil
              :initarg :reference
              :accessor reference-of
              :documentation "A reference sequence within which the
              interval lies.")
   (lower :type fixnum
          :initform 0
          :initarg :lower
          :accessor lower-of
          :documentation "The lower bound of the interval. If a
          reference is defined, this must be within the bounds of the
          reference sequence.")
   (upper :type fixnum
          :initform 0
          :initarg :upper
          :accessor upper-of
          :documentation "The upper bound of the interval."))
  (:documentation "An interval within a reference sequence. If a
  reference is defined, this must be within the bounds of the
  reference sequence. The basic interval has no notion of sequence
  strandedness; the bounds always refer to the forward strand."))

(defclass bio-sequence-interval (bio-sequence interval)
  ()
  (:documentation "A biological sequence that is an interval within a
  reference sequence."))

(defclass na-sequence-interval (na-sequence interval stranded-mixin)
  ()
  (:documentation "A nucleic acid sequence that is an interval within
  a reference sequence. In addition to the upper and lower bounds, a
  strand is defined. The strand indicates the strand of the reference
  sequence on which the interval lies."))

(defclass aa-sequence-interval (aa-sequence interval)
  ()
  (:documentation "An amino acid sequence that is an interval within
  a reference sequence."))

(defgeneric invert-complement (na-sequence-interval)
  (:documentation "Returns a copy of NA-SEQUENCE-INTERVAL at the
  corresponding position on the complementary strand of the reference
  sequence."))

(defgeneric ninvert-complement (na-sequence-interval)
  (:documentation "Destructively modifies NA-SEQUENCE-INTERVAL, moving
  it to the position on the complementary strand of the reference
  sequence."))

;;; Initialization methods
(defmethod initialize-instance :after ((interval bio-sequence-interval) &key)
  (with-slots (lower upper reference)
      interval
    (check-interval-range lower upper reference)))
  
(defmethod initialize-instance :after
    ((interval na-sequence-interval) &key)
  (with-slots (reference lower upper strand num-strands)
      interval
    (check-interval-range lower upper reference)
    (check-interval-strands strand num-strands reference)))

;;; Printing methods
(defmethod print-object ((interval interval) stream)
  (with-slots (lower upper)
      interval
    (format stream "#<INTERVAL ~a ~a>" lower upper)))

(defmethod print-object ((interval bio-sequence-interval) stream)
  (with-slots (lower upper)
      interval
    (format stream "#<BIO-SEQUENCE-INTERVAL ~a ~a>" lower upper)))

(defmethod print-object ((interval na-sequence-interval) stream)
  (with-slots (lower upper strand)
      interval
    (format stream "#<NA-SEQUENCE-INTERVAL ~a ~a ~a>" lower upper strand)))

;;; Implementation methods
(defmethod (setf reference-of) :before (value (interval interval))
  (with-slots (lower upper)
      interval
    (check-interval-range lower upper value)))

(defmethod (setf upper-of) :before (value (interval interval))
  (with-slots (reference lower)
      interval
    (check-interval-range lower value reference)))

(defmethod (setf lower-of) :before (value (interval interval))
  (with-slots (reference upper)
      interval
    (check-interval-range value upper reference)))

(defmethod (setf reference-of) :before
    (value (interval na-sequence-interval))
  (with-slots (strand num-strands)
      interval
    (check-interval-strands strand num-strands value)))

(defmethod (setf strand-of) :before
    (value (interval na-sequence-interval))
   (with-slots (reference num-strands)
       interval
     (check-interval-strands value num-strands reference)))

(defmethod (setf num-strands-of) :before
    (value (interval na-sequence-interval))
  (with-slots (reference strand)
      interval
    (check-interval-strands strand value reference)))

(defmethod length-of ((interval interval))
  (with-slots (lower upper)
      interval
    (- upper lower)))

(defmethod to-string :around ((interval interval) &key
                              start end token-case)
  (declare (ignore start end token-case))
  (with-slots (reference lower upper)
      interval
    (if reference
        (call-next-method)
      (make-string (- upper lower) :element-type 'base-char
                   :initial-element +gap-char+))))

(defmethod to-string ((interval interval) &key
                      (start 0) (end (length-of interval)) token-case)
  (with-slots (reference lower upper)
      interval
    (to-string reference :start (+ lower start) :end (+ lower end)
               :token-case token-case)))

(defmethod to-string ((interval na-sequence-interval) &key
                      (start 0) (end (length-of interval)) token-case)
  (with-slots (reference lower upper strand)
      interval
    (with-slots (num-strands)
        reference
      (cond ((or (eql *unknown-strand* strand)
                 (eql *forward-strand* strand))
             (to-string reference :start (+ lower start) :end (+ lower end)
                        :token-case token-case))
            ((and (eql *reverse-strand* strand) (= 2 num-strands))
             (to-string
              (nreverse-complement
               (subsequence reference (+ lower start) (+ lower end)))))
            (t
             (error 'bio-sequence-op-error
                    :text "a reverse-strand interval may not be created on a single-stranded sequence."))))))

(defmethod subsequence ((interval bio-sequence-interval) (start fixnum)
                        &optional end)
  (let ((end (or end (length-of interval))))
    (with-slots (reference lower upper)
        interval
      (subsequence reference (+ lower start) (+ lower end)))))

(defmethod subsequence ((interval na-sequence-interval)
                        (start fixnum) &optional end)
  (let ((end (or end (length-of interval))))
    (with-slots (reference lower upper strand num-strands)
        interval
      (let ((ref-num-strands (num-strands-of reference)))
        (cond ((or (eql *unknown-strand* strand)
                   (eql *forward-strand* strand))
               (subsequence reference (+ lower start) (+ lower end)))
              ((and (eql *reverse-strand* strand) (= ref-num-strands 2))
               (nreverse-complement
                (subsequence reference (+ lower start) (+ lower end))))
              (t
               (error 'bio-sequence-op-error
                      :text "a reverse-strand interval may not be applied to a single-stranded sequence.")))))))

(defmethod invert-complement ((interval na-sequence-interval))
  (with-slots (reference lower upper strand)
      interval
    (let ((ref-length (length-of reference))
          (int-length (- upper lower)))
      (make-instance 'na-sequence-interval
                     :lower (- ref-length upper)
                     :upper (+ (- ref-length upper) int-length)
                     :strand (invert-strand strand)
                     :reference reference))))

(defmethod ninvert-complement ((interval na-sequence-interval))
  (with-slots (reference lower upper strand num-strands)
      interval
    (check-interval-strands (invert-strand strand) num-strands reference)
    (let ((ref-length (length-of reference))
          (int-length (- upper lower)))
      (setf lower (- ref-length upper)
            upper (+ lower int-length)
            strand (invert-strand strand))))
  interval)

;;; Utility functions
(declaim (inline check-interval-range))
(defun check-interval-range (lower upper reference)
  "Validates interval bounds LOWER and UPPER against REFERENCE to
ensure that the interval lies within the bounds of the reference."
  (cond (reference
         (let ((length (length-of reference)))
           (cond ((or (< lower 0) (> lower length))
                  (error 'invalid-argument-error
                         :params 'lower
                         :args lower
                         :text "lower index must be >0 and <= sequence length"))
                 ((or (< upper 0) (> upper length))
                  (error 'invalid-argument-error
                         :params 'upper
                         :args upper
                         :text "upper index must be >0 and <= sequence length"))
                 ((< upper lower)
                  (error 'invalid-argument-error
                         :params '(lower upper)
                         :args (list lower upper)
                         :text "lower bound must be equal to or less than upper bound"))
                 (t t))))
          (t
           (cond ((< upper lower)
                  (error 'invalid-argument-error
                         :params '(lower upper)
                         :args (list lower upper)
                         :text "lower bound must be equal to or less than upper bound"))
                 (t t)))))

(declaim (inline check-interval-strands))
(defun check-interval-strands (strand num-strands reference)
  "Validates STRAND and NUM-STRANDS against REFERENCE to ensure that a
double-stranded interval is not applied to a single-stranded reference
and a reverse-strand interval is not applied to a single-stranded
reference."
  (if reference
      (let ((num-ref-strands (num-strands-of reference)))
        (cond ((and (= 1 num-ref-strands) (= 2 num-strands))
               (error 'invalid-argument-error
                      :params 'reference
                      :args reference
                      :text "a double-stranded interval may not be applied to a single-stranded sequence"))
              ((and (= 1 num-ref-strands) (eql *reverse-strand* strand))
                 (error 'invalid-argument-error
                        :params 'reference
                        :args reference
                        :text "a reverse-strand interval may not be applied to a single-stranded sequence"))
              (t t))))
  t)
