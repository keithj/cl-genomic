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

(defmacro define-subst-index (name elements &key (test '=))
  "Defines a numeric indexing function of ELEMENTS which accepts a
single argument. If the argument is equal \(by TEST\) to the nth
member of ELEMENTS, the function returns n. If the argument is not a
member of ELEMENTS, an error is raised."
  `(progn
    (defun ,name (value)
      (declare (optimize (speed 3) (safety 0)))
      (declare (type fixnum value))
      (cond ,@(loop
                 for elt in elements
                 for i = 0 then (1+ i)
                 collect `((,test value ,elt)
                           ,i))
            (t
             (error 'invalid-argument-error
                    :params 'value
                    :args value
                    :text "unknown subsitution matrix value"))))))

(defun combined-error-prob (e1 e2)
  "Calculates the combined error probability of a substitution between
two bases having individual error probabilities of E1 and E2, given
that the probability of a match occuring between two erronously called
bases is 1/3 * e1 * e2."
  (- (+ e1 e2) (* 4/3 (* e1 e2))))

(defun quality-match-score (e)
  "Calculates the match score given an error probabilty E."
  (let ((x (- 1.0 e)))
    (log (/ (if (zerop x)
                single-float-epsilon
              x) 0.25) 2)))

(defun quality-mismatch-penalty (e)
  "Calculates the mismatch penalty given an error probabilty E."
  (log (/ e 0.75) 2))

(defconstant +quality-match-matrix+
  (make-array '(100 100)
              :element-type 'single-float
              :initial-contents
              (loop
                 for q1 from 0 to 99
                 collect (loop
                            for q2 from 0 to 99
                            collect (quality-match-score
                                     (combined-error-prob
                                      (phred-probability q1)
                                      (phred-probability q2))))))
  "A sequence quality-adaptive DNA substitution matrix for matching
bases calculated as described by Malde in Bioinformatics 24,
pp. 897-900.")

(defconstant +quality-mismatch-matrix+
  (make-array '(100 100)
              :element-type 'single-float
              :initial-contents
              (loop
                 for q1 from 0 to 99
                 collect (loop
                            for q2 from 0 to 99
                            collect (quality-mismatch-penalty
                                     (combined-error-prob
                                      (phred-probability q1)
                                      (phred-probability q2))))))
  "A sequence quality-adaptive DNA substitution matrix for mismatching
bases calculated as described by Malde in Bioinformatics 24,
pp. 897-900.")

(defconstant +simple-dna-matrix+
  (make-array '(5 5)
              :element-type 'single-float
              :initial-contents '(( 5.0 -4.0 -4.0 -4.0 -1.0)
                                  (-4.0  5.0 -4.0 -4.0 -1.0)
                                  (-4.0 -4.0  5.0 -4.0 -1.0)
                                  (-4.0 -4.0 -4.0  5.0 -1.0)
                                  (-1.0 -1.0 -1.0 -1.0 -1.0)))
  "A DNA substitution matrix for sequences containing the bases A, C,
G, T and N.")

(defconstant +iupac-dna-matrix+
  (make-array '(15 15)
              :element-type 'single-float
              :initial-contents
              '(( 5.0 -4.0 -4.0 -4.0  2.0 -1.0  2.0  2.0 -1.0 -1.0  1.0  1.0  1.0 -2.0 -1.0)
                (-4.0  5.0 -4.0 -4.0 -1.0  2.0  2.0 -1.0  2.0 -1.0 -2.0  1.0  1.0  1.0 -1.0)
                (-4.0 -4.0  5.0 -4.0  2.0 -1.0 -1.0 -1.0  2.0  2.0  1.0 -2.0  1.0  1.0 -1.0)
                (-4.0 -4.0 -4.0  5.0 -1.0  2.0 -1.0  2.0 -1.0  2.0  1.0  1.0 -2.0  1.0 -1.0)
                ( 2.0 -1.0  2.0 -1.0  2.0 -2.0 -1.0  1.0  1.0  1.0  1.0 -1.0  1.0 -1.0 -1.0)
                (-1.0  2.0 -1.0  2.0 -2.0  2.0 -1.0  1.0  1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0)
                ( 2.0  2.0 -1.0 -1.0 -1.0 -1.0  2.0  1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0 -1.0)
                ( 2.0 -1.0 -1.0  2.0  1.0  1.0  1.0  2.0 -1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0)
                (-1.0  2.0  2.0 -1.0  1.0  1.0  1.0 -1.0  2.0  1.0 -1.0 -1.0  1.0  1.0 -1.0)
                (-1.0 -1.0  2.0  2.0  1.0  1.0 -1.0  1.0  1.0  2.0  1.0 -1.0 -1.0  1.0 -1.0)
                ( 1.0 -2.0  1.0  1.0  1.0 -1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0)
                ( 1.0  1.0 -2.0  1.0 -1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0  1.0 -1.0 -1.0 -1.0)
                ( 1.0  1.0  1.0 -2.0  1.0 -1.0  1.0 -1.0  1.0 -1.0 -1.0 -1.0  1.0 -1.0 -1.0)
                (-2.0  1.0  1.0  1.0 -1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0 -1.0 -1.0  1.0 -1.0)
                (-1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0)))
  "A DNA substitution matrix for sequences containing IUPAC bases.")

(defconstant +blosum50-matrix+
  (make-array '(23 23)
              :element-type 'single-float
              :initial-contents
              '(( 5.0 -2.0 -1.0 -2.0 -1.0 -1.0 -1.0  0.0 -2.0 -1.0 -2.0 -1.0 -1.0 -3.0 -1.0  1.0  0.0 -3.0 -2.0  0.0 -2.0 -1.0 -1.0)
                (-2.0  7.0 -1.0 -2.0 -4.0  1.0  0.0 -3.0  0.0 -4.0 -3.0  3.0 -2.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -1.0  0.0 -1.0)
                (-1.0 -1.0  7.0  2.0 -2.0  0.0  0.0  0.0  1.0 -3.0 -4.0  0.0 -2.0 -4.0 -2.0  1.0  0.0 -4.0 -2.0 -3.0  4.0  0.0 -1.0)
                (-2.0 -2.0  2.0  8.0 -4.0  0.0  2.0 -1.0 -1.0 -4.0 -4.0 -1.0 -4.0 -5.0 -1.0  0.0 -1.0 -5.0 -3.0 -4.0  5.0  1.0 -1.0)
                (-1.0 -4.0 -2.0 -4.0 13.0 -3.0 -3.0 -3.0 -3.0 -2.0 -2.0 -3.0 -2.0 -2.0 -4.0 -1.0 -1.0 -5.0 -3.0 -1.0 -3.0 -3.0 -2.0)
                (-1.0  1.0  0.0  0.0 -3.0  7.0  2.0 -2.0  1.0 -3.0 -2.0  2.0  0.0 -4.0 -1.0  0.0 -1.0 -1.0 -1.0 -3.0  0.0  4.0 -1.0)
                (-1.0  0.0  0.0  2.0 -3.0  2.0  6.0 -3.0  0.0 -4.0 -3.0  1.0 -2.0 -3.0 -1.0 -1.0 -1.0 -3.0 -2.0 -3.0  1.0  5.0 -1.0)
                ( 0.0 -3.0  0.0 -1.0 -3.0 -2.0 -3.0  8.0 -2.0 -4.0 -4.0 -2.0 -3.0 -4.0 -2.0  0.0 -2.0 -3.0 -3.0 -4.0 -1.0 -2.0 -2.0)
                (-2.0  0.0  1.0 -1.0 -3.0  1.0  0.0 -2.0 10.0 -4.0 -3.0  0.0 -1.0 -1.0 -2.0 -1.0 -2.0 -3.0  2.0 -4.0  0.0  0.0 -1.0)
                (-1.0 -4.0 -3.0 -4.0 -2.0 -3.0 -4.0 -4.0 -4.0  5.0  2.0 -3.0  2.0  0.0 -3.0 -3.0 -1.0 -3.0 -1.0  4.0 -4.0 -3.0 -1.0)
                (-2.0 -3.0 -4.0 -4.0 -2.0 -2.0 -3.0 -4.0 -3.0  2.0  5.0 -3.0  3.0  1.0 -4.0 -3.0 -1.0 -2.0 -1.0  1.0 -4.0 -3.0 -1.0)
                (-1.0  3.0  0.0 -1.0 -3.0  2.0  1.0 -2.0  0.0 -3.0 -3.0  6.0 -2.0 -4.0 -1.0  0.0 -1.0 -3.0 -2.0 -3.0  0.0  1.0 -1.0)
                (-1.0 -2.0 -2.0 -4.0 -2.0  0.0 -2.0 -3.0 -1.0  2.0  3.0 -2.0  7.0  0.0 -3.0 -2.0 -1.0 -1.0  0.0  1.0 -3.0 -1.0 -1.0)
                (-3.0 -3.0 -4.0 -5.0 -2.0 -4.0 -3.0 -4.0 -1.0  0.0  1.0 -4.0  0.0  8.0 -4.0 -3.0 -2.0  1.0  4.0 -1.0 -4.0 -4.0 -2.0)
                (-1.0 -3.0 -2.0 -1.0 -4.0 -1.0 -1.0 -2.0 -2.0 -3.0 -4.0 -1.0 -3.0 -4.0 10.0 -1.0 -1.0 -4.0 -3.0 -3.0 -2.0 -1.0 -2.0)
                ( 1.0 -1.0  1.0  0.0 -1.0  0.0 -1.0  0.0 -1.0 -3.0 -3.0  0.0 -2.0 -3.0 -1.0  5.0  2.0 -4.0 -2.0 -2.0  0.0  0.0 -1.0)
                ( 0.0 -1.0  0.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  2.0  5.0 -3.0 -2.0  0.0  0.0 -1.0  0.0)
                (-3.0 -3.0 -4.0 -5.0 -5.0 -1.0 -3.0 -3.0 -3.0 -3.0 -2.0 -3.0 -1.0  1.0 -4.0 -4.0 -3.0 15.0  2.0 -3.0 -5.0 -2.0 -3.0)
                (-2.0 -1.0 -2.0 -3.0 -3.0 -1.0 -2.0 -3.0  2.0 -1.0 -1.0 -2.0  0.0  4.0 -3.0 -2.0 -2.0  2.0  8.0 -1.0 -3.0 -2.0 -1.0)
                ( 0.0 -3.0 -3.0 -4.0 -1.0 -3.0 -3.0 -4.0 -4.0  4.0  1.0 -3.0  1.0 -1.0 -3.0 -2.0  0.0 -3.0 -1.0  5.0 -4.0 -3.0 -1.0)
                (-2.0 -1.0  4.0  5.0 -3.0  0.0  1.0 -1.0  0.0 -4.0 -4.0  0.0 -3.0 -4.0 -2.0  0.0  0.0 -5.0 -3.0 -4.0  5.0  2.0 -1.0)
                (-1.0  0.0  0.0  1.0 -3.0  4.0  5.0 -2.0  0.0 -3.0 -3.0  1.0 -1.0 -4.0 -1.0  0.0 -1.0 -2.0 -2.0 -3.0  2.0  5.0 -1.0)
                (-1.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -1.0  0.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0))))


;; "ACGTN" 4bit
(declaim (inline simple-dna-index))
(define-subst-index simple-dna-index
    (4 2 8 1 15))

;; "ACGTRYMWSKDHVBN" 4bit
(declaim (inline iupac-dna-index))
(define-subst-index iupac-dna-index
    (4 2 8 1 12 3 6 5 10 9 13 7 14 11 15))

;; "ARNDCQEGHILKMFPSTWYVBZX" 7bit
(declaim (inline blosum-50-index))
(define-subst-index blosum-50-index
    (1 18 14 4 3 17 5 7 8 9 12 11 13 6 16 19 20 23 25 22 2 26 24))

(declaim (inline simple-dna-subst))
(defun simple-dna-subst (x y)
  "Returns a substitution score from a simple DNA matrix \(permitted
residues are A, C, G, T and N\) for 4bit encoded DNA residues X and
Y."
  (aref +simple-dna-matrix+ (simple-dna-index x) (simple-dna-index y)))

(declaim (inline iupac-dna-subst))
(defun iupac-dna-subst (x y)
  "Returns a substitution score from a IUPAC DNA matrix \(permitted
residues are A, C, G, T, R, Y, M, W, S, K, D, H, V, B and N\) for 4bit
encoded DNA residues X and Y."
  (aref +iupac-dna-matrix+ (iupac-dna-index x) (iupac-dna-index y)))

(declaim (inline quality-dna-subst))
(defun quality-dna-subst (x y qx qy)
  (if (= x y)
      (aref +quality-match-matrix+ qx qy)
    (aref +quality-mismatch-matrix+ qx qy)))

(declaim (inline blosum-50-subst))
(defun blosum-50-subst (x y)
  (aref +blosum50-matrix+ (blosum-50-index x) (blosum-50-index y)))


