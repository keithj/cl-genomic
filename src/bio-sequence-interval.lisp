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

;;; Allen interval algebra
;;;
;;; x before y, y after x
;;; ------  -----
;;;
;;; x meets y, y met-by x
;;; ------
;;;       -----
;;;
;;; x overlaps y, y overlaps x
;;; ------
;;;     -----
;;;
;;; x starts y, y started-by x
;;; ------
;;; -----
;;;
;;; x during y, y contains x
;;;  ---
;;; -----
;;;
;;; x finishes y, y finished-by x
;;;   ---
;;; -----
;;;
;;; x equals y
;;; -----
;;; -----

(defgeneric beforep (x y)
  (:documentation "Returns T if X is before Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric afterp (x y)
  (:documentation "Returns T if X is after Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric meetsp (x y)
  (:documentation "Returns T if X meets Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric met-by-p (x y)
  (:documentation "Returns T if X is met by Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric startsp (x y)
  (:documentation "Returns T if X starts Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric started-by-p (x y)
  (:documentation "Returns T if X is started by Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric duringp (x y)
  (:documentation "Returns T if X occurs during Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric containsp (x y)
  (:documentation "Returns T if X contains Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric finishesp (x y)
  (:documentation "Returns T if X finishes Y according to Allen's
  Interval Algebra, or NIL otherwise."))

(defgeneric finished-by-p (x y)
  (:documentation "Returns T if X is finished by Y according to
  Allen's Interval Algebra, or NIL otherwise."))

(defgeneric interval-equal (x y)
  (:documentation "Returns T if X is interval equal Y according to
  Allen's Interval Algebra, or NIL otherwise."))

;;; Initialization methods
(defmethod initialize-instance :after ((interval na-sequence-interval) &key)
  (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of)
                   (strand strand-of) (num-strands num-strands-of))
      interval
    (%check-interval-range lower upper reference)
    (%check-interval-strands strand num-strands reference)))

(defmethod initialize-instance :after ((interval aa-sequence-interval) &key)
  (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of))
      interval
    (%check-interval-range lower upper reference)))

;;; Printing methods
(defmethod print-object ((interval interval) stream)
  (with-accessors ((lower lower-of) (upper upper-of))
      interval
    (format stream "#<INTERVAL ~a ~a>" lower upper)))

(defmethod print-object ((interval na-sequence-interval) stream)
  (with-accessors ((lower lower-of) (upper upper-of) (strand strand-of))
      interval
    (format stream "#<NA-SEQUENCE-INTERVAL ~a ~a ~a>" lower upper strand)))

(defmethod print-object ((interval aa-sequence-interval) stream)
  (with-accessors ((lower lower-of) (upper upper-of))
      interval
    (format stream "#<AA-SEQUENCE-INTERVAL ~a ~a>" lower upper)))

;;; Implementation methods
(defmethod (setf reference-of) :before (value (interval interval))
  (with-accessors ((lower lower-of) (upper upper-of))
      interval
    (%check-interval-range lower upper value)))

(defmethod (setf upper-of) :before (value (interval interval))
  (with-accessors ((lower lower-of) (reference reference-of))
      interval
    (%check-interval-range lower value reference)))

(defmethod (setf lower-of) :before (value (interval interval))
  (with-accessors ((upper upper-of) (reference reference-of))
      interval
    (%check-interval-range value upper reference)))

(defmethod (setf reference-of) :before
    (value (interval na-sequence-interval))
  (with-accessors ((strand strand-of) (num-strands num-strands-of))
      interval
    (%check-interval-strands strand num-strands value)))

(defmethod (setf strand-of) :before
    (value (interval na-sequence-interval))
  (with-accessors ((reference reference-of) (num-strands num-strands-of))
       interval
    (%check-interval-strands value num-strands reference)))

(defmethod (setf num-strands-of) :before
    (value (interval na-sequence-interval))
  (with-accessors ((reference reference-of) (strand strand-of))
      interval
    (%check-interval-strands strand value reference)))

(defmethod simplep ((interval na-sequence-interval))
  (with-accessors ((reference reference-of))
      interval
    (and reference (simplep reference))))

(defmethod simplep ((interval aa-sequence-interval))
  (with-accessors ((reference reference-of))
      interval
    (and reference (simplep reference))))

(defmethod ambiguousp ((interval na-sequence-interval))
  (with-accessors ((reference reference-of))
      interval
    (or (null reference) (ambiguousp reference))))

(defmethod ambiguousp ((interval aa-sequence-interval))
  (with-accessors ((reference reference-of))
      interval
    (or (null reference) (ambiguousp reference))))

(defmethod virtualp ((interval na-sequence-interval))
  (with-accessors ((reference reference-of))
      interval
    (or (null reference) (virtualp reference))))

(defmethod virtualp ((interval aa-sequence-interval))
  (with-accessors ((reference reference-of))
      interval
    (or (null reference) (virtualp reference))))

(defmethod length-of ((interval interval))
  (with-accessors ((lower lower-of) (upper upper-of))
      interval
    (- upper lower)))

(defmethod coerce-sequence :around ((interval interval) (type (eql 'string))
                                    &key start end)
  (declare (ignore start end))
  (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of))
      interval
    (if reference
        (call-next-method)
      (make-string (- upper lower) :element-type 'base-char
                   :initial-element +gap-char+))))

(defmethod coerce-sequence ((interval interval) (type (eql 'string))
                            &key (start 0) (end (length-of interval)))
  (with-accessors ((reference reference-of))
      interval
    (coerce-sequence reference 'string :start start :end end)))

(defmethod coerce-sequence ((interval na-sequence-interval)
                            (type (eql 'string))
                            &key (start 0) (end (length-of interval)))
  (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of)
                   (strand strand-of))
      interval
    (with-accessors ((num-strands num-strands-of))
        reference
      (cond ((or (eql *unknown-strand* strand)
                 (eql *forward-strand* strand))
             (coerce-sequence reference 'string :start (+ lower start)
                              :end (+ lower end)))
            ((and (eql *reverse-strand* strand) (= 2 num-strands))
             (coerce-sequence
              (nreverse-complement
               (subsequence reference (+ lower start) (+ lower end))) 'string))
            (t
             (error 'bio-sequence-op-error
                    :text "a reverse-strand interval may not be created on a single-stranded sequence."))))))

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

(defmethod subsequence ((interval aa-sequence-interval) (start fixnum)
                        &optional end)
  (let ((end (or end (length-of interval))))
    (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of))
        interval
      (subsequence reference (+ lower start) (+ lower end)))))

(defmethod invert-complement ((interval na-sequence-interval))
  (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of)
                   (strand strand-of))
      interval
    (let ((ref-length (length-of reference))
          (int-length (- upper lower)))
      (make-instance 'na-sequence-interval
                     :lower (- ref-length upper)
                     :upper (+ (- ref-length upper) int-length)
                     :strand (invert-strand strand)
                     :reference reference))))

(defmethod ninvert-complement ((interval na-sequence-interval))
  (with-accessors ((reference reference-of)
                   (strand strand-of) (num-strands num-strands-of))
      interval
    (%check-interval-strands (invert-strand strand) num-strands reference)
    (with-slots (lower upper) ; direct slot access to allow inversion in place
        interval
      (let ((ref-length (length-of reference))
            (int-length (- upper lower)))
        (setf lower (- ref-length upper)
              upper (+ lower int-length)
              strand (invert-strand strand)))))
  interval)

(defmethod beforep ((x interval) (y interval))
  (< (upper-of x) (lower-of y)))

(defmethod afterp ((x interval) (y interval))
  (> (lower-of x) (upper-of y)))

(defmethod meetsp ((x interval) (y interval))
  (= (upper-of x) (lower-of y)))

(defmethod met-by-p ((x interval) (y interval))
  (= (lower-of x) (upper-of y)))

(defmethod startsp ((x interval) (y interval))
  (and (= (lower-of x) (lower-of y))
       (< (upper-of x) (upper-of y))))

(defmethod started-by-p ((x interval) (y interval))
  (and (= (lower-of x) (lower-of y))
       (> (upper-of x) (upper-of y))))

(defmethod duringp ((x interval) (y interval))
  (and (>= (lower-of x) (lower-of y))
       (<= (upper-of x) (upper-of y))))

(defmethod containsp ((x interval) (y interval))
  (and (<= (lower-of x) (lower-of y))
       (>= (upper-of x) (upper-of y))))

(defmethod finishesp ((x interval) (y interval))
  (and (> (lower-of x) (lower-of y))
       (= (upper-of x) (upper-of y))))

(defmethod finished-by-p ((x interval) (y interval))
  (and (< (lower-of x) (lower-of y))
       (= (upper-of x) (upper-of y))))

(defmethod interval-equal ((x interval) (y interval))
  (and (= (lower-of x) (lower-of y))
       (= (upper-of x) (upper-of y))))


;; FIXME -- circular sequences

;; These may have intervals that travel around the sequence multiple
;; times

;; TODO --
;; canonicalize negative coordinates to positive (maybe public function?)
;; calculate numbers of rotations for intervals (maybe public function?)
;; overlaps, intersections, unions


;;; Utility functions
(declaim (inline %check-interval-range))
(defun %check-interval-range (lower upper reference)
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

(declaim (inline %check-interval-strands))
(defun %check-interval-strands (strand num-strands reference)
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

