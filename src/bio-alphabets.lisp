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
   (tokens :initform ""
           :initarg :tokens
           :reader tokens-of
           :documentation "The set of member tokens of the
alphabet.")
   (index :initarg :index
          :reader index-of
          :documentation "An index of the residues in the alphabet."))
  (:documentation "Alphabets are sets of tokens."))

(defconstant +gap-char+ #\-)
(defconstant +encoded-gap-char+ #b0000)

(defvar *dna*
  (let ((tokens (make-array 16
                            :element-type 'base-char
                            :initial-contents "acgtrykmswbdhvn-")))
    (make-instance 'alphabet
                   :name :dna
                   :tokens tokens
                   :index
                   (loop
                      with index = (make-hash-table)
                      for i from 0 below (length tokens)
                      do (setf (gethash
                                (encode-dna-4bit (aref tokens i)) index) i)
                      finally (return index))))
  "The IUPAC DNA alphabet.")

(defvar *rna*
  (let ((tokens (make-array 16
                            :element-type 'base-char
                            :initial-contents "acgurykmswbdhvn-")))
    (make-instance 'alphabet
                   :name :rna
                   :tokens tokens
                   :index
                   (loop
                      with index = (make-hash-table)
                      for i from 0 below (length tokens)
                      do (setf (gethash
                                (encode-rna-4bit (aref tokens i)) index) i)
                      finally (return index))))
  "The IUPAC RNA alphabet.")

(defvar *aa*
  (let ((tokens (make-array 28
                            :element-type 'base-char
                            :initial-contents "ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")))
    (make-instance 'alphabet
                   :name :aa
                   :tokens tokens
                   :index
                   (loop
                      with index = (make-hash-table)
                      for i from 0 below (length tokens)
                      do (setf (gethash
                                (encode-aa-7bit (aref tokens i)) index) i)
                      finally (return index)))))

(defvar *alphabets* (make-hash-table)
  "The standard biological alphabets.")

(setf (gethash :dna *alphabets*) *dna*)
(setf (gethash :rna *alphabets*) *rna*)
(setf (gethash :aa *alphabets*) *aa*)
