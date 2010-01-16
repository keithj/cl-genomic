;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
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

(defclass sequence-strand ()
  ((name :initarg :name
         :reader name-of
         :documentation "The strand name.")
   (symbol-of :initarg :symbol
              :reader symbol-of
              :documentation "The Lisp symbol representing the strand.")
   (token :initarg :token
          :reader token-of
          :documentation "The token representing the strand.")
   (number :initarg :number
           :reader number-of
           :documentation "The number representing the strand."))
  (:documentation "The strand of a nucleotide sequence."))

(defvar *forward-strand*
  (make-instance 'sequence-strand :name "forward" :symbol :forward
                 :token #\+ :number 1)
  "The forward nucleotide strand.")

(defvar *reverse-strand*
  (make-instance 'sequence-strand :name "reverse" :symbol :reverse
                 :token #\- :number -1)
  "The reverse nucleotide strand.")

(defvar *unknown-strand*
  (make-instance 'sequence-strand :name "unknown" :symbol :unknown
                 :token #\? :number nil)
  "An unknown nucleotide strand.")

(defclass stranded-mixin ()
  ((strand :type sequence-strand
           :initform *unknown-strand*
           :initarg :strand
           :reader strand-of
           :documentation "The nucleotide strand.")))

(defclass identity-mixin ()
  ((identity :initform nil
             :initarg :identity
             :reader identity-of
             :documentation "A temporary, locally unique identifier."))
  (:documentation "A mixin which allows assignment of a temporary
local identifier to an object. An identity of NIL signifies an
anonymous object."))

(defclass description-mixin ()
  ((description :initform nil
                :initarg :description
                :accessor description-of
                :documentation "A free text description."))
  (:documentation "A mixin which allows assignment of a free text
description to an object."))

(defclass quality-mixin ()
  ((metric :initform (error "A metric is required.")
           :initarg :metric
           :reader metric-of
           :documentation "A description of the quality metric
measured by the quality values. For example, p-value, Phred score or
Illumina score. This should be changed to a controlled vocabulary or
enumeration.")
   (quality :initform (error "Quality values are required.")
            :initarg :quality
            :accessor quality-of
            :documentation "The array of quality values which should
be the same length as the array of residues."))
  (:documentation "A mixin with support for bio-sequences that have a
numeric quality value for each residue."))

(defclass token-sequence ()
  ((alphabet :initarg :alphabet
             :reader alphabet-of
             :documentation "The alphabet of the sequence."))
  (:documentation "A sequence of tokens belonging to an alphabet."))

(defclass virtual-token-sequence (token-sequence)
  ((length :initform 0
           :initarg :length
           :accessor length-of
           :documentation "The length of the sequence."))
  (:documentation "A token sequence whose tokens are not explicitly
represented."))

(defclass encoded-token-sequence (token-sequence)
  ()
  (:documentation "A token sequence whose tokens have been transformed
to a representation other than characters, typically small integers."))

(defclass encoded-vector-sequence (encoded-token-sequence)
  ((vector :initform (make-array 0 :element-type 'fixnum)
           :initarg :vector
           :accessor vector-of
           :documentation "The token vector of the sequence."))
  (:documentation "An encoded sequence whose tokens are stored in an
in-memory vector."))

(defclass bio-sequence ()
  ()
  (:documentation "A biological sequence."))

(defclass na-sequence (bio-sequence)
  ((num-strands :initform 2
                :initarg :num-strands
                :accessor num-strands-of
                :documentation "The number of sequence strands."))
  (:documentation "A nucleic acid sequence."))

(defclass dna-sequence (na-sequence)
  ((alphabet :allocation :class
             :initform *dna*))
  (:documentation "A concrete DNA sequence comprising IUPAC ambiguity
bases."))

(defclass rna-sequence (na-sequence)
  ((alphabet :allocation :class
             :initform *rna*))
  (:documentation "A concrete RNA sequence comprising IUPAC ambiguity
bases."))

(defclass aa-sequence (bio-sequence)
  ((alphabet :allocation :class
             :initform *aa*)))

(defclass encoded-dna-sequence (dna-sequence encoded-vector-sequence
                                identity-mixin description-mixin)
  ())

(defclass encoded-rna-sequence (rna-sequence encoded-vector-sequence
                                identity-mixin description-mixin)
  ())

(defclass dna-quality-sequence (encoded-dna-sequence quality-mixin
                                identity-mixin description-mixin)
  ())

(defclass virtual-dna-sequence (dna-sequence virtual-token-sequence
                                identity-mixin description-mixin)
  ())

(defclass virtual-rna-sequence (rna-sequence virtual-token-sequence
                                identity-mixin description-mixin)
  ())

(defclass encoded-aa-sequence (aa-sequence encoded-vector-sequence
                               identity-mixin description-mixin)
  ())

(defclass virtual-aa-sequence (aa-sequence virtual-token-sequence
                               identity-mixin description-mixin)
  ())

(defclass mapped-dna-sequence (dna-sequence token-sequence
                               dxn:mapped-vector-char
                               identity-mixin description-mixin)
  ())

;; Also need to add circularity
