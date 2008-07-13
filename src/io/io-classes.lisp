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
  ()
  (:documentation "The base class of all biological sequence
parsers. The default methods specialised on this class are all
no-ops, ignoring any data and returning NIL."))

(defclass quality-parser-mixin ()
  ((metric :initform nil
           :initarg :metric
           :reader parsed-metric
           :documentation "The quality metric, e.g. :phred, :illumina,
parsed from the input stream."))
  (:documentation "A parser specialised for processing biological
sequence data with additional residue quality information."))

(defclass raw-sequence-parser (bio-sequence-parser
                               quality-parser-mixin)
  ((raw :initform '()
        :accessor parsed-raw
        :documentation "The raw sequence data parsed from the input
stream."))
  (:documentation "A parser specialised for processing raw biological
sequence data. This class is typically used for simple reformatting,
splitting or counting operations where making CLOS objects is not
desirable."))

(defclass simple-sequence-parser (bio-sequence-parser)
  ((alphabet :initform nil
             :accessor parsed-alphabet
             :documentation "The sequence alphabet designator,
e.g. :dna, :rna, parsed from the input stream.")
   (identity :initform nil
             :accessor parsed-identity
             :documentation "The sequence identity parsed from the
input stream.")
   (description :initform nil
                :accessor parsed-description
                :documentation "The sequence documentation parsed from
the input stream.")
   (residues :initform (make-array 0 :adjustable t :fill-pointer 0)
             :accessor parsed-residues
             :documentation "The sequence residues parsed from the
input stream."))
  (:documentation "A parser specialised for processing biological
sequence data to build CLOS objects."))

(defclass quality-sequence-parser (simple-sequence-parser
                                   quality-parser-mixin)
  ((quality :initform (make-array 0 :adjustable t :fill-pointer 0)
            :reader parsed-quality
            :documentation "The sequence quality data parsed from the
input stream."))
  (:documentation "A parser specialised for processing biological
sequence data with quality to build CLOS objects."))

(defclass virtual-sequence-parser (simple-sequence-parser)
  ((length :initform 0
           :accessor parsed-length
           :documentation "The sequence length parsed from the input
stream."))
  (:documentation "A parser specialised for processing biological
sequence data to build CLOS objects that do not contain explicit
residue data."))
