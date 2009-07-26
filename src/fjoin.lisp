;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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
                                        :upper most-positive-fixnum))

(defun sentinelp (interval)
  (eql *sentinel* interval))

(defun make-interval-gen (array)
  (let* ((i 0)
         (j (length array))
         (current (aref array i)))
    (defgenerator
        :current current
        :next (prog1
                  current
                (setf i (1+ i)
                      current (if (< i j)
                                  (aref array i)
                                *sentinel*)))
        :more t)))

(defun fjoin (ix iy test)
  (let ((wx (make-queue))
        (wy (make-queue))
        (overlaps (make-array 100 :adjustable t :fill-pointer 0)))
    (loop
       with x = (next ix) 
       with y = (next iy)
       while (not (and (sentinelp x) (sentinelp y)))
       do (cond ((<= (lower-of x) (lower-of y))
                 (dolist (overlap (scan test x wx y wy))
                   (vector-push-extend overlap overlaps))
                 (setf x (next ix)))
                (t
                 (dolist (overlap (scan test y wy x wx))
                   (vector-push-extend overlap overlaps))
                 (setf y (next iy))))
       finally (return overlaps))))

(defun scan (test f wf g wg)
  (when (not (inclusive-beforep f g))
    (queue-enqueue f wf))
  (loop
     for h in (queue-head (queue-dequeue-if (lambda (x)
                                              (inclusive-beforep x f)) wg))
     when (funcall test h f)            ; (inclusive-overlapsp h f)
     collect h))
