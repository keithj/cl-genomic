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

(defconstant +gap-char+ #\-
  "The gap character.")
(defconstant +encoded-gap-char+ #b0000
  "The encoded gap character in all encoded alphabets.")
(defconstant +terminator-char+ #\*
  "The peptide sequence terminator character.")

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
          :documentation "An index of the tokens in the alphabet."))
  (:documentation "Alphabets are sets of tokens."))

(defun make-alphabet (name tokens encoder)
  (make-instance 'alphabet
                 :name name
                 :tokens tokens
                 :index (loop
                           with index = (make-hash-table)
                           for i from 0 below (length tokens)
                           do (setf (gethash
                                     (funcall encoder (elt tokens i)) index) i)
                           finally (return index))))

(defun make-nth-order-alphabet (name alphabet n encoder)
  (let ((tokens (permutations-of-n (tokens-of alphabet) n)))
    (make-instance 'alphabet
                   :name name
                   :tokens tokens
                   :index (loop
                             with index = (make-hash-table :test #'equal)
                             for i from 0 below (length tokens)
                             do (setf (gethash
                                       (mapcar encoder (elt tokens i)) index) i)
                             finally (return index)))))

(defun permutations-of-n (elements &optional (n 1))
  (let ((m (1- n)))
    (if (zerop m)
        (mapcar #'list elements)
      (let ((perms ()))
        (dolist (e elements)
          (dolist (p (permutations-of-n elements m))
            (push (cons e p) perms)))
        (nreverse perms)))))

(defvar *simple-dna*
  (let ((tokens '(#\a #\c #\g #\t)))
    (make-alphabet :simple-dna tokens #'encode-dna-4bit))
  "The simple DNA alphabet.")

(defvar *simple-rna*
  (let ((tokens '(#\a #\c #\g #\u)))
    (make-alphabet :simple-dna tokens #'encode-rna-4bit))
  "The simple RNA alphabet.")

(defvar *dna*
  (let ((tokens `(#\a #\c #\g #\t #\r #\y #\k #\m #\s #\w
                  #\b #\d #\h #\v #\n ,+gap-char+)))
    (make-alphabet :dna tokens #'encode-dna-4bit))
  "The IUPAC DNA alphabet.")

(defvar *rna*
  (let ((tokens `(#\a #\c #\g #\u #\r #\y #\k #\m #\s #\w
                  #\b #\d #\h #\v #\n ,+gap-char+)))
    (make-alphabet :rna tokens #'encode-rna-4bit))
  "The IUPAC RNA alphabet.")

(defvar *aa*
  (let ((tokens `(#\A #\B #\C #\D #\E #\F #\G #\H #\I #\J #\K
                  #\L #\M #\N #\O #\P #\Q #\R #\S #\T #\U #\V
                  #\W #\X #\Y #\Z ,+terminator-char+ ,+gap-char+)))
    (make-alphabet :aa tokens #'encode-aa-7bit))
  "The amino acid alphabet.")

(defvar *dna-codons*
  (make-nth-order-alphabet :dna-codons *simple-dna* 3 #'encode-dna-4bit)
  "The alphabet of 64 standard codons.")

(defvar *rna-codons*
  (make-nth-order-alphabet :rna-codons *simple-rna* 3 #'encode-rna-4bit)
  "The alphabet of 64 standard codons.")

(defvar *alphabets* (make-hash-table)
  "The standard biological alphabets.")

(setf (gethash :simple-dna *alphabets*) *simple-dna*)
(setf (gethash :simple-rna *alphabets*) *simple-rna*)
(setf (gethash :dna *alphabets*) *dna*)
(setf (gethash :rna *alphabets*) *rna*)
(setf (gethash :aa *alphabets*) *aa*)
(setf (gethash :dna-codons *alphabets*) *dna-codons*)
(setf (gethash :rna-codons *alphabets*) *rna-codons*)
