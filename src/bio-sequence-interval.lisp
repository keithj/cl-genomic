;;;
;;; Copyright (C) 2008-2009 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
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
          :reader lower-of
          :documentation "The lower bound of the interval. If a
          reference is defined, this must be within the bounds of the
          reference sequence.")
   (upper :type fixnum
          :initform 0
          :initarg :upper
          :reader upper-of
          :documentation "The upper bound of the interval."))
  (:documentation "An interval within a reference sequence. If a
reference is defined, this must be within the bounds of the reference
sequence. The basic interval has no notion of sequence strandedness;
the bounds always refer to the forward strand. Spatial relationships
between intervals are described using Allen interval algebra.

 Allen interval algebra (x is the upper interval)

   x before y, y after x
   ------
          -----

   x meets y, y met-by x (neither x or y are zero-width)
   ------
         -----

   x overlaps y, y overlaps x
   ------
       -----

   x starts y, y started-by x (neither x or y are zero-width)
   ----
   -------

   x during y, y contains x
    ---
   -----

   x finishes y, y finished-by x (neither x or y are zero-width)
     ---
   -----

   x equals y
   -----
   -----

This set of functions is extended with others that are less strict.

  inclusive-before which is the union of before and meets
  inclusive-after which is the union of after and met-by
  inclusive-contains which is the union of contains, started-by, finished-by
    and equals
  inclusive-overlaps which is the union of overlaps and inclusive-contains"))

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

(defgeneric beforep (x y)
  (:documentation "Returns T if X is before Y according to Allen's
Interval Algebra, or NIL otherwise. This definition of 'before' is
stricter than is often used in bioinformatics. See also
{defun inclusive-beforep} ."))

(defgeneric afterp (x y)
  (:documentation "Returns T if X is after Y according to Allen's
Interval Algebra, or NIL otherwise. This definition of 'after' is
stricter than is often used in bioinformatics. See also
{defun inclusive-afterp} ."))

(defgeneric meetsp (x y)
  (:documentation "Returns T if X meets Y according to Allen's
Interval Algebra, or NIL otherwise."))

(defgeneric met-by-p (x y)
  (:documentation "Returns T if X is met by Y according to Allen's
Interval Algebra, or NIL otherwise."))

(defgeneric overlapsp (x y)
  (:documentation "Returns T if X overlaps Y according to Allen's
Interval Algebra, or NIL otherwise. This definition of 'overlaps' is
stricter than is often used in bioinformatics. See also
{defun inclusive-overlapsp} ."))

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
Interval Algebra, or NIL otherwise. This definition of 'contains' is
stricter than is often used in bioinformatics. See also
{defun inclusive-containsp} ."))

(defgeneric finishesp (x y)
  (:documentation "Returns T if X finishes Y according to Allen's
Interval Algebra, or NIL otherwise."))

(defgeneric finished-by-p (x y)
  (:documentation "Returns T if X is finished by Y according to
Allen's Interval Algebra, or NIL otherwise."))

(defgeneric interval-equal (x y)
  (:documentation "Returns T if X is interval equal Y according to
Allen's Interval Algebra, or NIL otherwise."))

(defgeneric inclusive-beforep (x y)
  (:documentation "The union of beforep and meetsp. This predicate is
often named 'before' in bioinformatics use cases. See also
{defun beforep} ."))

(defgeneric inclusive-afterp (x y)
  (:documentation "The union of afterp and meet-by-p. This predicate is
often named 'after' in bioinformatics use cases. See also
{defun beforep} ."))

(defgeneric inclusive-containsp (x y)
  (:documentation "The union of containsp, started-by-p, finished-by-p
and interval-equal. This predicate is often named 'contains' in
bioinformatics use cases. See also {defun containsp} ."))

(defgeneric inclusive-duringp (x y)
  (:documentation "The union of duringp, startsp, finishesp and
interval-equal. See also {defun duringp} ."))

(defgeneric inclusive-overlapsp (x y)
  (:documentation "The union of overlapsp, inclusive-duringp and
inclusive-containsp. This predicate is often named 'overlaps' in
bioinformatics use cases. See also {defun overlapsp} ."))

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

(defmethod make-interval ((reference na-sequence) &rest initargs)
  (apply #'make-instance 'na-sequence-interval :reference reference initargs))

(defmethod make-interval ((reference aa-sequence) &rest initargs)
  (apply #'make-instance 'aa-sequence-interval :reference reference initargs))

;; (let* ((ref (make-dna (make-string 10000000 :initial-element #\a)))
;;        (dna (find-alphabet :simple-dna))
;;        (intervals (make-array 100000 :initial-element nil)))
;;   (with-sequence-residues (residue ref)
;;     (setf residue (random-token-of dna)))
;;   (dotimes (i (length intervals))
;;     (let* ((len (* 10 (random 2000)))
;;            (lower (random (- 10000000 len))))
;;       (setf (svref intervals i)
;;             (make-interval ref :lower lower :upper (+ lower len)))))
;;   (sort intervals #'< :key #'lower-of)
;;   (time (let ((x (svref intervals 0)))
;;           (cons x
;;                 (loop
;;                    for i across intervals
;;                    when (overlaps x i)
;;                    collect i)))))

;;; Implementation methods
(defmethod (setf reference-of) :before (value (interval interval))
  (with-accessors ((lower lower-of) (upper upper-of))
      interval
    (%check-interval-range lower upper value)))

(defmethod (setf reference-of) :before (value (interval na-sequence-interval))
  (with-accessors ((strand strand-of) (num-strands num-strands-of))
      interval
    (%check-interval-strands strand num-strands value)))

(defmethod (setf num-strands-of) :before (value (interval na-sequence-interval))
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
      (cond ((or (eql *unknown-strand* strand) (eql *forward-strand* strand))
             (coerce-sequence reference 'string :start (+ lower start)
                              :end (+ lower end)))
            ((and (eql *reverse-strand* strand) (= 2 num-strands))
             (coerce-sequence
              (nreverse-complement
               (subsequence reference (+ lower start) (+ lower end))) 'string))
            (t
             (error 'bio-sequence-op-error
                    :text (txt "a reverse-strand interval may not be created"
                               "on a single-stranded sequence.")))))))

(defmethod subsequence ((interval na-sequence-interval)
                        (start fixnum) &optional end)
  (let ((end (or end (length-of interval))))
    (with-slots (reference lower upper strand num-strands)
        interval
      (let ((ref-num-strands (num-strands-of reference)))
        (cond ((or (eql *unknown-strand* strand) (eql *forward-strand* strand))
               (subsequence reference (+ lower start) (+ lower end)))
              ((and (eql *reverse-strand* strand) (= ref-num-strands 2))
               (nreverse-complement
                (subsequence reference (+ lower start) (+ lower end))))
              (t
               (error 'bio-sequence-op-error
                      :text (txt "a reverse-strand interval may not be applied"
                                 "to a single-stranded sequence."))))))))

(defmethod subsequence ((interval aa-sequence-interval) (start fixnum)
                        &optional end)
  (let ((end (or end (length-of interval))))
    (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of))
        interval
      (subsequence reference (+ lower start) (+ lower end)))))

(defmethod reverse-complement ((interval na-sequence-interval))
  (with-accessors ((lower lower-of) (upper upper-of) (reference reference-of)
                   (strand strand-of))
      interval
    (let ((ref-length (length-of reference))
          (int-length (- upper lower)))
      (make-instance 'na-sequence-interval
                     :lower (- ref-length upper)
                     :upper (+ (- ref-length upper) int-length)
                     :strand (complement-strand strand)
                     :reference reference))))

(defmethod nreverse-complement ((interval na-sequence-interval))
  (with-accessors ((reference reference-of)
                   (strand strand-of) (num-strands num-strands-of))
      interval
    (%check-interval-strands (complement-strand strand) num-strands reference)
    (with-slots (lower upper strand) ; direct slot access to allow
                                     ; inversion in place
        interval
      (let ((ref-length (length-of reference))
            (int-length (- upper lower)))
        (setf lower (- ref-length upper)
              upper (+ lower int-length)
              strand (complement-strand strand)))))
  interval)

(defmethod beforep ((x interval) (y interval))
  (< (upper-of x) (lower-of y)))

(defmethod afterp ((x interval) (y interval))
  (> (lower-of x) (upper-of y)))

(defmethod meetsp ((x interval) (y interval))
  (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (and (< lowx upx)             ; to meet y, x may not be zero-width
           (= upx lowy)))))

(defmethod met-by-p ((x interval) (y interval))
  (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (and (< lowy upy)             ; to meet x, y may not be zero-width
           (= lowx upy)))))

(defmethod overlapsp ((x interval) (y interval))
  (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (or (< lowx lowy upx upy)
          (< lowy lowx upy upx)))))

(defmethod startsp ((x interval) (y interval))
   (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (and (< lowx upx) (< lowy upy)   ; neither x or y are zero-width
           (= lowx lowy) (< upx upy)))))

(defmethod started-by-p ((x interval) (y interval))
  (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (and (< lowx upx) (< lowy upy)   ; neither x or y are zero-width
           (= lowx lowy) (> upx upy)))))

(defmethod duringp ((x interval) (y interval))
  (and (> (lower-of x) (lower-of y))
       (< (upper-of x) (upper-of y))))

(defmethod containsp ((x interval) (y interval))
  (and (< (lower-of x) (lower-of y))
       (> (upper-of x) (upper-of y))))

(defmethod finishesp ((x interval) (y interval))
  (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (and (< lowx upx) (< lowy upy)   ; neither x or y are zero-width
           (> lowx lowy) (= upx upy)))))

(defmethod finished-by-p ((x interval) (y interval))
  (with-accessors ((lowx lower-of) (upx upper-of))
      x
    (with-accessors ((lowy lower-of) (upy upper-of))
        y
      (and (< lowx upx) (< lowy upy)   ; neither x or y are zero-width
           (< lowx lowy) (= upx upy)))))

(defmethod interval-equal ((x interval) (y interval))
  (and (= (lower-of x) (lower-of y))
       (= (upper-of x) (upper-of y))))

(defmethod inclusive-beforep ((x interval) (y interval))
  (or (beforep x y) (meetsp x y)))

(defmethod inclusive-afterp ((x interval) (y interval))
  (or (afterp x y) (met-by-p x y)))

(defmethod inclusive-containsp ((x interval) (y interval))
  (or (interval-equal x y) (containsp x y) 
      (started-by-p x y) (finished-by-p x y)))

(defmethod inclusive-duringp ((x interval) (y interval))
  (or (interval-equal x y) (duringp x y)
      (startsp x y) (finishesp x y)))

(defmethod inclusive-overlapsp ((x interval) (y interval))
  (or (overlapsp x y) (inclusive-duringp x y) (inclusive-containsp x y)))


;; TODO -- separate methods or add an integer parameter to existing
;; overlapsp method?

;; (defmethod min-overlapsp ((x interval) (y interval) (n fixnum))
;;   (with-accessors ((lowx lower-of) (upx upper-of))
;;       x
;;     (with-accessors ((lowy lower-of) (upy upper-of))
;;         y
;;       (>= (- (min lowx lowy) (max upx upy)) n))))

;; (defmethod min-beforep ((x interval) (y interval) (n fixnum))
;;   (<= (upper-of x) (+ (lower-of y) n)))


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
                         :params '(lower length)
                         :args (list lower length)
                         :text "lower must be >0 and <= length"))
                 ((or (< upper 0) (> upper length))
                  (error 'invalid-argument-error
                         :params '(upper length)
                         :args (list upper length)
                         :text "upper must be >0 and <= length"))
                 ((< upper lower)
                  (error 'invalid-argument-error
                         :params '(lower upper)
                         :args (list lower upper)
                         :text "lower must be <= upper"))
                 (t t))))
          (t
           (cond ((< upper lower)
                  (error 'invalid-argument-error
                         :params '(lower upper)
                         :args (list lower upper)
                         :text "lower must be <= upper"))
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
                      :text (txt "a double-stranded interval may not be"
                                 "applied to a single-stranded sequence")))
              ((and (= 1 num-ref-strands) (eql *reverse-strand* strand))
                 (error 'invalid-argument-error
                        :params 'reference
                        :args reference
                        :text (txt "a reverse-strand interval may not be"
                                   "applied to a single-stranded sequence")))
              (t t))))
  t)
