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

(defvar *alphabets* (make-hash-table))

(defun find-alphabet (name)
  (multiple-value-bind (alphabet presentp)
      (gethash name *alphabets*)
    (unless presentp
      (error "Invalid alphabet ~a." name))
    alphabet))

(defclass alphabet ()
  ((name :initarg :name
         :reader name-of
         :documentation "The alphabet name.")
   (encoder :initarg :encoder
            :reader encoder-of)
   (decoder :initarg :decoder
            :reader decoder-of)
   (encoded-index :initarg :encoded-index
                  :reader encoded-index-of)
   (tokens :initform ""
             :initarg :tokens
             :reader tokens-of
             :documentation "The set of member tokens of the
alphabet."))
  (:documentation "Alphabets are sets of tokens."))

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

(setf (gethash :dna *alphabets*)
      (make-instance 'alphabet
                     :name :dna
                     :encoder #'encode-dna-2bit
                     :decoder #'decode-dna-2bit
                     :tokens (make-array 4
                                         :element-type 'base-char
                                         :initial-contents "acgt")))
(setf (gethash :rna *alphabets*)
      (make-instance 'alphabet
                     :name :rna
                     :encoder #'encode-rna-2bit
                     :decoder #'decode-rna-2bit
                     :tokens (make-array 4
                                         :element-type 'base-char
                                         :initial-contents "acgu")))
(setf (gethash :iupac-dna *alphabets*)
      (make-instance 'alphabet
                     :name :iupac-dna
                     :encoder #'encode-dna-4bit
                     :decoder #'decode-dna-4bit
                     :tokens
                     (make-array 15
                                 :element-type 'base-char
                                 :initial-contents "acgtrykmswbdhvn")))
(setf (gethash :iupac-rna *alphabets*)
      (make-instance 'alphabet
                     :name :iupac-rna
                     :encoder #'encode-rna-4bit
                     :decoder #'decode-rna-4bit
                     :tokens
                     (make-array 15
                                 :element-type 'base-char
                                 :initial-contents "acgurykmswbdhvn")))

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
local identifier to an object."))

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
be the same length as the array of residue tokens."))
  (:documentation "A mixin with support for bio-sequences that have a
numeric quality value for each residue."))

(defclass bio-sequence (identity-mixin)
  ((alphabet :initarg :alphabet
             :accessor alphabet-of
             :documentation "The alphabet whose tokens comprise the
sequence.")
   (token-seq :initform nil
              :initarg :token-seq
              :documentation "The residue tokens of the sequence or
NIL if no sequence data are available.")
   (length :initform nil
           :initarg :length
           :documentation "The logical length of the sequence in
residues. This value is not required to be supported by concrete
sequence data."))
  (:documentation "A biological sequence comprising tokens from a
specified alphabet. Its position relative to other bio-sequences may
be given by specifying the sequence-locations of its start in those
sequences."))

(defclass nucleic-acid-sequence (bio-sequence)
  ()
  (:documentation "A logical nucleic acid sequence."))

(defclass dna-sequence (nucleic-acid-sequence)
  ()
  (:documentation "A logical DNA sequence."))

(defclass rna-sequence (nucleic-acid-sequence)
  ()
  (:documentation "A logical RNA sequence."))

(defclass iupac-dna-sequence (dna-sequence)
  ((alphabet :initform (find-alphabet :iupac-dna)
             :allocation :class))
  (:documentation "A concrete DNA sequence comprising IUPAC ambiguity
bases."))

(defclass simple-dna-sequence (dna-sequence)
  ((alphabet :initform (find-alphabet :dna)
             :allocation :class))
  (:documentation "A concrete DNA sequence comprising unambiguous
bases T, C, A and G."))

(defclass iupac-rna-sequence (rna-sequence)
  ((alphabet :initform (find-alphabet :iupac-rna)
             :allocation :class))
  (:documentation "A concrete RNA sequence comprising IUPAC ambiguity
bases."))

(defclass simple-rna-sequence (rna-sequence)
  ((alphabet :initform (find-alphabet :rna)
             :allocation :class))
  (:documentation "A concrete RNA sequence comprising unambiguous
bases U, C, A and G."))

(defclass simple-dna-quality-sequence (simple-dna-sequence quality-mixin)
  ())

(defclass iupac-dna-quality-sequence (iupac-dna-sequence quality-mixin)
  ())
