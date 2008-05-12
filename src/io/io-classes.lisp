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

(defclass bio-sequence-parser ()
  ())

(defclass quality-parser-mixin ()
  ((metric :initform nil
           :initarg :metric
           :reader parsed-metric)))

(defclass raw-sequence-parser (bio-sequence-parser
                               quality-parser-mixin)
  ((raw :initform '()
        :accessor parsed-raw)))

(defclass simple-sequence-parser (bio-sequence-parser)
  ((alphabet :initform nil
             :accessor parsed-alphabet)
   (identity :initform nil
             :accessor parsed-identity)
   (description :initform nil
                :accessor parsed-description)
   (residues :initform (make-array 0 :adjustable t :fill-pointer 0)
             :accessor parsed-residues)))

(defclass quality-sequence-parser (simple-sequence-parser
                                   quality-parser-mixin)
  ((quality :initform (make-array 0 :adjustable t :fill-pointer 0)
            :reader parsed-quality)))

(defclass virtual-sequence-parser (simple-sequence-parser)
  ((length :initform 0
           :accessor parsed-length)))
