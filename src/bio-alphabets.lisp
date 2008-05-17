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
   (decoded-index :initarg :decoded-index
                  :reader decoded-index-of)
   (tokens :initform ""
           :initarg :tokens
           :reader tokens-of
           :documentation "The set of member tokens of the
alphabet."))
  (:documentation "Alphabets are sets of tokens."))

(defvar *dna*
  (make-instance 'alphabet
                 :name :dna
                 :encoder #'encode-dna-4bit
                 :decoder #'decode-dna-4bit
                 :tokens (make-array 15
                                     :element-type 'base-char
                                     :initial-contents "acgtrykmswbdhvn"))
  "The IUPAC DNA alphabet.")
(defvar *rna*
  (make-instance 'alphabet
                 :name :rna
                 :encoder #'encode-rna-4bit
                 :decoder #'decode-rna-4bit
                 :tokens (make-array 15
                                     :element-type 'base-char
                                     :initial-contents "acgurykmswbdhvn"))
  "The IUPAC RNA alphabet.")

(defvar *alphabets* (make-hash-table)
  "The standard biological alphabets.")

(setf (gethash :dna *alphabets*) *dna*)
(setf (gethash :rna *alphabets*) *rna*)
