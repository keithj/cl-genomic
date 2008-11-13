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

(defstruct nucleotide-base
  (name "" :type string)
  (description "" :type string)
  (token #\n :type base-char))

(defconstant +gap-char+ #\-
  "The gap character.")
(defconstant +encoded-gap-char+ #b0000
  "The encoded gap character in all encoded alphabets.")
(defconstant +codon-size+ 3
  "The number of bases in a codon.")

(defvar *alphabets* (make-hash-table)
  "The standard biological alphabets.")

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
  "Returns a new simple alphabet where the tokens are atomic. Use this
function for creating alphabets such as DNA and RNA. Use
{defun make-nth-order-alphabet} to create alphabets of compound tokens,
such as the codon alphabet.

Arguments:

- name (symbol): The alphabet symbolic name.
- tokens (list atoms): The alphabet tokens.
- encoder (function): A function that accepts a single alphabet token
  and returns its encoded value.

Returns:

- An {defclass alphabet} ."
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
  "Returns a new alphabet where the tokens are lists of tokens from a
simple alphabet. Use this function for creating alphabets of compound
tokens, such as the codon alphabet.

Arguments:

- name (symbol): The alphabet symbolic name.
- alphabet ( {defclass alphabet} ): The alphabet whose tokens will
  be grouped to form the tokens of the new alphabet.
- n (integer): The number of tokens of the simple alphabet that
  combine to form a token of the new alphabet.
- encoder (function): A function that accepts a single simple alphabet
  token and returns its encoded value.

Returns:

- An {defclass alphabet} ."
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

(defun register-alphabet (alphabet)
  "Registers a global standard ALPHABET."
  (setf (gethash (name-of alphabet) *alphabets*) alphabet))

(defun find-alphabet (name)
  "Returns a standard ALPHABET designated by a symbol NAME, such
as :dna :rna or :aa."
  (multiple-value-bind (alphabet presentp)
      (gethash name *alphabets*)
    (unless presentp
      (error 'invalid-argument-error
             :params 'name
             :args name
             :text "no such alphabet"))
    alphabet))

(defun registered-alphabets ()
  "Returns a list of all registered alphabets."
  (loop
     for alphabet being the hash-values of *alphabets*
     collect alphabet into alphabets
     finally (return (sort alphabets #'string<=
                           :key (lambda (alpha)
                                  (symbol-name (name-of alpha)))))))

(defun permutations-of-n (elements &optional (n 1))
  "Returns a list of permutations of N ELEMENTS."
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
                  #\W #\X #\Y #\Z #\* ,+gap-char+)))
    (make-alphabet :aa tokens #'encode-aa-7bit))
  "The amino acid alphabet.")

(defvar *dna-codons*
  (make-nth-order-alphabet :dna-codons *simple-dna* 3 #'encode-dna-4bit)
  "The alphabet of 64 standard codons.")

(defvar *rna-codons*
  (make-nth-order-alphabet :rna-codons *simple-rna* 3 #'encode-rna-4bit)
  "The alphabet of 64 standard codons.")

(register-alphabet *simple-dna*)
(register-alphabet *simple-rna*)
(register-alphabet *dna*)
(register-alphabet *rna*)
(register-alphabet *aa*)
(register-alphabet *dna-codons*)
(register-alphabet *rna-codons*)
