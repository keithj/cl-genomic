;;;
;;; Copyright (C) 2009-2010 Keith James. All rights reserved.
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

(defparameter *sentinel* (make-instance 'na-sequence-interval
                                        :lower most-positive-fixnum
                                        :upper most-positive-fixnum)
  "The sentinel interval used by the fjoin algorithm.")

(defun sentinelp (interval)
  "Returns T if INTERVAL is the sentinel, or NIL otherwise."
  (eql *sentinel* interval))

(defmethod make-interval-input ((vector vector))
  "Returns a new generator that will return the elements of VECTOR,
followed by *SENTINEL*."
  (let* ((i 0)
         (j (length vector))
         (current (aref vector i)))
    (defgenerator
        :current current
        :next (prog1
                  current
                (setf i (1+ i)
                      current (if (< i j)
                                  (aref vector i)
                                *sentinel*)))
        :more t)))

(defmacro with-interval-input ((var obj) &body body)
  `(let ((,var (make-interval-input ,obj)))
    ,@body))

(defgeneric compare-intervals (query-intervals reference-intervals
                               exclude test callback)
  (:documentation "Compares QUERY-INTERVALS with REFERENCE-INTERVALS
using he pairwise fjoin alogorithm described in
Richardson, J. E. (2006) Journal of Computational Biology 13 (8),
1457-1464.

Arguments:

- query-intervals (sequence intervals): A sequence of intervals.
- reference-intervals (sequence intervals): A sequence of intervals.
- exclude (function): A function of two arguments used to test
  intervals for eviction from the working window. In the fjoin paper
  this is leftOf. The cl-genomic equivalent is {defun beforep} .
- test (function): A function of two arguments used to test intervals
  to determine whether they are related. This function should return T
  if the realtionship holds, or NIL otherwise. In the fjoin paper this
  is overlaps. The cl-genomic equivalent is {defun inclusive-overlapsp} .
- callback (function): A function of two arguments that is called
  whenever a pair of intervals are found to relate by function
  TEST. CALLBACK is called with the two intervals as arguments.

Returns:

- T."))

(defmethod compare-intervals ((x vector) (y vector) exclude test callback)
  (with-interval-input (ix x)
    (with-interval-input (iy y)
      (fjoin ix iy exclude test callback))))

(defun fjoin (ix iy exclude test callback)
  "This function implements the pairwise fjoin alogorithm described in
Richardson, J. E. (2006) Journal of Computational Biology 13 (8),
1457-1464.

All arguments are functions in order to allow extension beyond simple
overlap detection. For example, detecting degree of overlap in terms
of percentage length or number of bases, degree of separation or
strand relationships.

Extension from pairwise relationships to a full n-ary fjoin is
planned.

This implementation is full of side-effects because it's a very
literal translation of the pseudocode in the fjoin paper. This is
intentional.

Arguments:

- ix (function): An iterator function that returns
 {defclass bio-sequence-interval} s in the order of their position on
a reference sequence, until exhausted when it returns the sentinel
interval.
- iy (function): An iterator function that returns
 {defclass bio-sequence-interval} s in the order of their position on
a reference sequence, until exhausted when it returns the sentinel
interval.
- exclude (function): A function of two arguments used to test
  intervals for eviction from the working window. In the fjoin paper
  this is leftOf. The cl-genomic equivalent is {defun beforep} .
- test (function): A function of two arguments used to test intervals
  to determine whether they are related. This function should return T
  if the realtionship holds, or NIL otherwise. In the fjoin paper this
  is overlaps. The cl-genomic equivalent is {defun inclusive-overlapsp} .
- callback (function): A function of two arguments that is called
  whenever a pair of intervals are found to relate by function
  TEST. CALLBACK is called with the two intervals as arguments.

Returns:

- T."
  (let ((wx (make-queue))
        (wy (make-queue)))
    (loop
       with x = (next ix)
       with y = (next iy)
       while (not (and (sentinelp x) (sentinelp y)))
       do (cond ((<= (lower-of x) (lower-of y))
                 (scan x wx y wy exclude test callback)
                 (setf x (next ix)))
                (t
                 (scan y wy x wx exclude test callback)
                 (setf y (next iy)))))
    t))

(defun scan (f wf g wg exclude test &optional callback)
  "Window-scanning part of the fjoin algorithm.

Arguments:

- f (interval): The current interval in the lagging stream.
- wf (queue): The window across the lagging stream.
- g (interval): The current interval in the leading stream.
- wf (queue): The window across the leading stream.
- exclude (function): A function of two arguments used to test
  intervals for eviction from the working window.
- test (function): A function of two arguments used to test intervals
  to determine whether they are related. This function should return T
  if the realtionship holds, or NIL otherwise.

Optional:

- callback (function): A function of two arguments that is called
  whenever a pair of intervals are found to relate by function
  TEST. CALLBACK is called with the two intervals as arguments."
  (unless (funcall exclude f g)
    (queue-enqueue f wf))
  (loop
     for x in (queue-head (queue-dequeue-if (lambda (y)
                                              (funcall exclude y f)) wg))
     do (when (and callback (funcall test x f))
          (funcall callback x f))))
