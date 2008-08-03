;;;
;;; Copyright (C) 2008 Keith James. All rights reserved.
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

(defvar *simple-dna-matrix*
  (make-array '(5 5)
              :element-type 'single-float
              :initial-contents '(( 5.0 -4.0 -4.0 -4.0 -1.0)
                                  (-4.0  5.0 -4.0 -4.0 -1.0)
                                  (-4.0 -4.0  5.0 -4.0 -1.0)
                                  (-4.0 -4.0 -4.0  5.0 -1.0)
                                  (-1.0 -1.0 -1.0 -1.0 -1.0))))

(defvar *simple-dna-index*
  (make-array 5 :element-type '(unsigned-byte 4)
              :initial-contents
              (loop
                 for c across "ACGTN"
                 collect (encode-dna-4bit c))))

(defvar *simple-dna-subst*
  (make-instance 'subst-matrix
                 :matrix *simple-dna-matrix*
                 :index *simple-dna-index*))

(defvar *blosum50-matrix*
  (make-array '(23 23)
              :element-type 'single-float
              :initial-contents
              '(( 5.0 -2.0 -1.0 -2.0 -1.0 -1.0 -1.0  0.0 -2.0 -1.0 -2.0 -1.0 -1.0 -3.0 -1.0  1.0  0.0 -3.0 -2.0  0.0 -2.0 -1.0 -1.0)
                (-2.0  7.0 -1.0 -2.0 -4.0  1.0  0.0 -3.0  0.0 -4.0 -3.0  3.0 -2.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -1.0  0.0 -1.0)
                (-1.0 -1.0  7.0  2.0 -2.0  0.0  0.0  0.0  1.0 -3.0 -4.0  0.0 -2.0 -4.0 -2.0  1.0  0.0 -4.0 -2.0 -3.0  4.0  0.0 -1.0)
                (-2.0 -2.0  2.0  8.0 -4.0  0.0  2.0 -1.0 -1.0 -4.0 -4.0 -1.0 -4.0 -5.0 -1.0  0.0 -1.0 -5.0 -3.0 -4.0  5.0  1.0 -1.0)
                (-1.0 -4.0 -2.0 -4.0 13.0 -3.0 -3.0 -3.0 -3.0 -2.0 -2.0 -3.0 -2.0 -2.0 -4.0 -1.0 -1.0 -5.0 -3.0 -1.0 -3.0 -3.0 -2.0)
                (-1.0  1.0  0.0  0.0 -3.0  7.0  2.0 -2.0  1.0 -3.0 -2.0  2.0  0.0 -4.0 -1.0  0.0 -1.0 -1.0 -1.0 -3.0  0.0  4.0 -1.0)
                (-1.0  0.0  0.0  2.0 -3.0  2.0  6.0 -3.0  0.0 -4.0 -3.0  1.0 -2.0 -3.0 -1.0 -1.0 -1.0 -3.0 -2.0 -3.0  1.0  5.0 -1.0)
                ( 0.0 -3.0  0.0 -1.0 -3.0 -2.0 -3.0  8.0 -2.0 -4.0 -4.0 -2.0 -3.0 -4.0 -2.0  0.0 -2.0 -3.0 -3.0 -4.0 -1.0 -2.0 -2.0)
                (-2.0  0.0  1.0 -1.0 -3.0  1.0  0.0 -2.0 10.0 -4.0 -3.0  0.0 -1.0 -1.0 -2.0 -1.0 -2.0 -3.0  2.0 -4.0  0.0  0.0 -1.0)
                (-1.0 -4.0 -3.0 -4.0 -2.0 -3.0 -4.0 -4.0 -4.0  5.0  2.0 -3.0  2.0  0.0 -3.0 -3.0 -1.0 -3.0 -1.0  4.0 -4.0 -3.0 -1.0)
                (-2.0 -3.0 -4.0 -4.0 -2.0 -2.0 -3.0 -4.0 -3.0  2.0  5.0 -3.0  3.0  1.0 -4.0 -3.0 -1.0 -2.0 -1.0  1.0 -4.0 -3.0 -1.0)
                (-1.0  3.0  0.0 -1.0 -3.0  2.0  1.0 -2.0  0.0 -3.0 -3.0  6.0 -2.0 -4.0 -1.0  0.0 -1.0 -3.0 -2.0 -3.0  0.0  1.0 -1.0)
                (-1.0 -2.0 -2.0 -4.0 -2.0  0.0 -2.0 -3.0 -1.0  2.0  3.0 -2.0  7.0  0.0 -3.0 -2.0 -1.0 -1.0  0.0  1.0 -3.0 -1.0 -1.0)
                (-3.0 -3.0 -4.0 -5.0 -2.0 -4.0 -3.0 -4.0 -1.0  0.0  1.0 -4.0  0.0  8.0 -4.0 -3.0 -2.0  1.0  4.0 -1.0 -4.0 -4.0 -2.0)
                (-1.0 -3.0 -2.0 -1.0 -4.0 -1.0 -1.0 -2.0 -2.0 -3.0 -4.0 -1.0 -3.0 -4.0 10.0 -1.0 -1.0 -4.0 -3.0 -3.0 -2.0 -1.0 -2.0)
                ( 1.0 -1.0  1.0  0.0 -1.0  0.0 -1.0  0.0 -1.0 -3.0 -3.0  0.0 -2.0 -3.0 -1.0  5.0  2.0 -4.0 -2.0 -2.0  0.0  0.0 -1.0)
                ( 0.0 -1.0  0.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  2.0  5.0 -3.0 -2.0  0.0  0.0 -1.0  0.0)
                (-3.0 -3.0 -4.0 -5.0 -5.0 -1.0 -3.0 -3.0 -3.0 -3.0 -2.0 -3.0 -1.0  1.0 -4.0 -4.0 -3.0 15.0  2.0 -3.0 -5.0 -2.0 -3.0)
                (-2.0 -1.0 -2.0 -3.0 -3.0 -1.0 -2.0 -3.0  2.0 -1.0 -1.0 -2.0  0.0  4.0 -3.0 -2.0 -2.0  2.0  8.0 -1.0 -3.0 -2.0 -1.0)
                ( 0.0 -3.0 -3.0 -4.0 -1.0 -3.0 -3.0 -4.0 -4.0  4.0  1.0 -3.0  1.0 -1.0 -3.0 -2.0  0.0 -3.0 -1.0  5.0 -4.0 -3.0 -1.0)
                (-2.0 -1.0  4.0  5.0 -3.0  0.0  1.0 -1.0  0.0 -4.0 -4.0  0.0 -3.0 -4.0 -2.0  0.0  0.0 -5.0 -3.0 -4.0  5.0  2.0 -1.0)
                (-1.0  0.0  0.0  1.0 -3.0  4.0  5.0 -2.0  0.0 -3.0 -3.0  1.0 -1.0 -4.0 -1.0  0.0 -1.0 -2.0 -2.0 -3.0  2.0  5.0 -1.0)
                (-1.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -1.0  0.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0))))

(defvar *blosum50-index*
  (make-array 23 :element-type '(unsigned-byte 8)
              :initial-contents "ARNDCQEGHILKMFPSTWYVBZX"))

(defvar *blosum50-subst*
  (make-instance 'subst-matrix
                 :matrix *blosum50-matrix*
                 :index *blosum50-index*))
