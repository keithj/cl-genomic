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

(in-package :cl-genomic-test)

(deftestsuite fjoin-tests (cl-genomic-tests)
  ())

(defun make-intervals (coord-pairs)
  (make-array (length coord-pairs)
              :initial-contents (mapcar
                                 (lambda (pair)
                                   (make-instance 'interval
                                                  :lower (first pair)
                                                  :upper (second pair)))
                                 coord-pairs)))

(defun test-fjoin-intervals (xspecs yspecs rspecs exclude test)
  (let ((x (make-intervals xspecs))
        (y (make-intervals yspecs))
        (expected-result (mapcar #'make-intervals rspecs))
        (result ()))
    (flet ((store-intervals (a b)
             (push (list a b) result)))
      (compare-intervals x y exclude test #'store-intervals))
    (setf result (reverse result))
    (ensure (every (lambda (observed expected)
                     (every #'interval-equal observed expected))
                   result expected-result)
            :report "Expected ~a but observed ~a"
            :arguments (expected-result result))))

;; TODO -- extend tests to cover all the interval predicates
(addtest (fjoin-tests) compare-intervals/1
  (let ((in1  '(( 0 9) (10 19) (20 29) (30 39) (40 49)))
        (in2  '(       (10 18) (20 28)         (40 48)))
        (out1 '(((10 19) (10 18))
                ((20 29) (20 28))
                ((40 49) (40 48))))
        (out2 '(((10 18) (10 19))
                ((20 28) (20 29))
                ((40 48) (40 49)))))
    (test-fjoin-intervals in1 in2 out1
                          #'beforep #'inclusive-overlapsp)
    (test-fjoin-intervals in2 in1 out2
                          #'beforep #'inclusive-overlapsp)))

(addtest (fjoin-tests) compare-intervals/2
  (let ((in1  '((0  9) (10 19) (20 29) (30 39) (40 49)))
        (in2  '((8         18)                 (48 59)))
        (out1 '(((0   9) ( 8 18))
                (( 8 18) (10 19))
                ((40 49) (48 59))))
        (out2 '(((0   9) ( 8 18))
                (( 8 18) (10 19))
                ((40 49) (48 59)))))
    (test-fjoin-intervals in1 in2 out1
                          #'beforep #'inclusive-overlapsp)
     (test-fjoin-intervals in2 in1 out2
                           #'beforep #'inclusive-overlapsp)))
