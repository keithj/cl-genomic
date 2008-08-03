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

(defclass sequence-strand ()
  ((name :initarg :name
         :reader name-of
         :documentation "The strand name.")
   (token :initarg :token
            :reader token-of
            :documentation "The token representing the strand.")
   (number :initarg :number
           :reader number-of
             :documentation "The number representing the strand."))
  (:documentation "The strand of a nucleotide sequence."))

(defvar *forward-strand*
  (make-instance 'sequence-strand
                 :name :forward
                 :token #\+
                 :number 1))

(defvar *reverse-strand*
  (make-instance 'sequence-strand
                 :name :reverse
                 :token #\-
                 :number -1))

(defvar *without-strand*
  (make-instance 'sequence-strand
                 :name :unstranded
                 :token #\.
                 :number 0))

(defvar *unknown-strand*
  (make-instance 'sequence-strand
                 :name :unknown
                 :token #\?
                 :number nil))

(defclass identity-mixin ()
  ((identity :initform nil
             :initarg :identity
             :accessor identity-of
             :documentation "A temporary locally unique identifier."))
  (:documentation "A mixin which allows assignment of a temporary
local identifier to an object. An identity of NIL signifies an
anonymous object."))

(defclass quality-mixin ()
  ((metric :initform (error "A metric is required.")
           :initarg :metric
           :reader metric-of
           :documentation "A description of the quality metric
measured by the quality values. For example, p-value, Phred score or
Illumina score. This should be changed to a controlled vocabulary or
enumeration.")
   (quality :initform (error "A quality argument is required.")
            :initarg :quality
            :accessor quality-of
            :documentation "The array of quality values which should
be the same length as the array of residues."))
  (:documentation "A mixin with support for bio-sequences that have a
numeric quality value for each residue."))

(defclass token-sequence ()
  ((alphabet :initarg :alphabet
             :reader alphabet-of))
  (:documentation "A sequence of tokens belonging to an alphabet."))

(defclass virtual-token-sequence (token-sequence)
  ((length :initarg :length
           :accessor length-of)))

(defclass simple-token-sequence (token-sequence)
  ())

(defclass encoded-token-sequence (token-sequence)
  ())

(defclass mmapped-sequence ()
  ())

(defclass vector-sequence ()
  ((vector :initarg :vector
           :accessor vector-of)))

(defclass simple-vector-sequence (simple-token-sequence vector-sequence)
  ())

(defclass encoded-vector-sequence (encoded-token-sequence vector-sequence)
  ())

(defclass bio-sequence ()
  ()
  (:documentation "A biological sequence."))

(defclass nucleic-acid-sequence (bio-sequence)
  ()
  (:documentation "A nucleic acid sequence."))

(defclass dna-sequence (nucleic-acid-sequence)
  ((alphabet :allocation :class
             :initform *dna*))
  (:documentation "A concrete DNA sequence comprising IUPAC ambiguity
bases."))

(defclass rna-sequence (nucleic-acid-sequence)
  ((alphabet :allocation :class
             :initform *rna*))
  (:documentation "A concrete RNA sequence comprising IUPAC ambiguity
bases."))

(defclass encoded-dna-sequence (dna-sequence encoded-vector-sequence
                                             identity-mixin)
  ())

(defclass encoded-rna-sequence (rna-sequence encoded-vector-sequence
                                             identity-mixin)
  ())

(defclass dna-quality-sequence (encoded-dna-sequence quality-mixin
                                                     identity-mixin)
  ())

(defclass virtual-dna-sequence (dna-sequence virtual-token-sequence
                                             identity-mixin)
  ())

(defclass virtual-rna-sequence (rna-sequence virtual-token-sequence
                                             identity-mixin)
  ())

(defclass subst-matrix ()
  ((matrix :initarg :matrix
           :reader matrix-of)
   (index :initarg :index
          :reader index-of)))
