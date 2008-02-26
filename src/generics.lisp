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

(defgeneric simplep (bio-sequence alphabet-designator)
  (:documentation "Returns T if BIO-SEQUENCE contains no ambiguous
residues according to the alphabets specified by ALPHABET-DESIGNATOR,
or NIL otherwise."))

(defgeneric virtualp (bio-sequence)
  (:documentation "Returns T if BIO-SEQUENCE has no concrete token-seq
representation, merely a length, or NIL otherwise."))

(defgeneric copy-sequence (bio-sequence)
  (:documentation "Returns a copy of BIO-SEQUENCE."))

(defgeneric length-of (bio-sequence)
  (:documentation "Returns the length of BIO-SEQUENCE."))

(defgeneric token-seq-of (bio-sequence)
  (:documentation "Returns the token sequence representing the
residues of BIO-SEQUENCE."))

(defgeneric residue-of (bio-sequence index)
  (:documentation "Returns the residue at INDEX of BIO-SEQUENCE."))

(defgeneric (setf residue-of) (value bio-sequence index)
  (:documentation "Sets the residue at INDEX of BIO-SEQUENCE to
VALUE."))

(defgeneric cumulative-lengths (ranges)
  (:documentation "Returns a list of integers which are the
cumulative total lengths of RANGES."))

(defgeneric to-string (bio-sequence &optional start length)
  (:documentation "Returns the string representing BIO-SEQUENCE,
starting at the first residue, or index START, and continuing to the
last residue, or for LENGTH residues."))


;;; bio-sequence io generics

(defgeneric read-bio-sequence (stream format &key alphabet ambiguity
                               virtualp &allow-other-keys)
  (:documentation "Reads a sequence record from STREAM of FORMAT
(e.g. :fasta, :fastq). Keywords are used to specify the expected
alphabet (:dna, :rna) and ambiguity (:iupac, nil). The VIRTUALP
keyword, if T, indicates that a virtual sequence should be created,
having the correct length, but no concrete residues. If no sequence
can be read, NIL should be returned. This is the high-level sequence
reading interface which returns bio-sequence CLOS objects."))

(defgeneric read-bio-sequence-alist (stream format &key alphabet ambiguity
                                     virtualp callback callback-args)
  (:documentation "Reads a sequence record of ALPHABET with AMBIGUITY
from STREAM in FORMAT, optionally applying function CALLBACK with
additional CALLBACK-ARGS to the result. The VIRTUALP keyword, if T,
indicates that a virtual sequence should be created, having the
correct length, but no concrete residues. This is the low-level
sequence reading interface which normally returns an alist.  If no
sequence can be read, NIL should be returned. An optional CALLBACK may
be supplied which should be a function accepting the alist as the
first argument, plus any number of additional arguments to be supplied
in the list CALLBACK-ARGS. The result of applying the callback should
be returned."))
