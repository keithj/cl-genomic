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


;;; bio-sequence generics

(defgeneric copy-sequence (bio-sequence)
  (:documentation "Returns a copy of BIO-SEQUENCE."))

(defgeneric length-of (bio-sequence)
  (:documentation "Returns the length of BIO-SEQUENCE."))

(defgeneric residue-of (bio-sequence index)
  (:documentation "Returns the residue at INDEX of BIO-SEQUENCE."))

(defgeneric (setf residue-of) (value bio-sequence index)
  (:documentation "Sets the residue at INDEX of BIO-SEQUENCE to
VALUE."))

(defgeneric to-string (bio-sequence &optional start length)
  (:documentation "Returns the string representing BIO-SEQUENCE,
starting at the first residue, or index START, and continuing to the
last residue, or for LENGTH residues."))


;;; bio-sequence io generics

(defgeneric read-bio-sequence (stream &key alphabet ambiguity format)
  (:documentation "Reads a sequence record from STREAM. Keywords are
used to specify the expected alphabet (:dna, :rna), ambiguity (:iupac,
nil) and record format (:fasta, :fastq)."))

(defgeneric read-bio-sequence-alist (input format alphabet ambiguity
                                     &optional callback callback-args)
  (:documentation "Reads a sequence record of ALPHABET with AMBIGUITY
from INPUT in FORMAT, optionally applying function CALLBACK with
additional CALLBACK-ARGS to the result."))
