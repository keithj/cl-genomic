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

(in-package :cl-genomic-test)

(deftestsuite bio-sequence-alignment-tests (cl-genomic-tests)
  ())

(defparameter *swg-aa-score*
  #2A((0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0)
      (0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0)
      (0.0  0.0  0.0  5.0  0.0  5.0  0.0  0.0  0.0  0.0  0.0)
      (0.0  0.0  0.0  0.0  2.0  0.0 20.0  0.0  0.0  0.0  0.0)
      (0.0 10.0  0.0  0.0  0.0  0.0  0.0 18.0 20.0  9.0  8.0)
      (0.0  0.0 16.0  0.0  0.0  0.0  0.0  7.0 18.0 26.0 16.0)
      (0.0  0.0  0.0 21.0  6.0 10.0  1.0  9.0  6.0 17.0 25.0)
      (0.0  0.0  6.0  5.0 18.0 10.0  7.0  6.0  9.0 15.0 23.0)))

(defparameter *swg-aa-insertx*
  #2A((0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0)
      (0.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0)
      (0.0 -10.0 -10.0 -10.0  -5.0  -6.0  -5.0  -6.0  -7.0  -8.0  -9.0)
      (0.0 -10.0 -10.0 -10.0 -10.0  -8.0  -9.0  10.0   9.0   8.0   7.0)
      (0.0 -10.0   0.0  -1.0  -2.0  -3.0  -4.0   0.0   8.0  10.0   9.0)
      (0.0 -10.0 -10.0   6.0   5.0   4.0   3.0   2.0   1.0   8.0  16.0)
      (0.0 -10.0 -10.0  -4.0  11.0  10.0   9.0   8.0   7.0   6.0   7.0)
      (0.0 -10.0 -10.0  -4.0   1.0   8.0   7.0   6.0   5.0   4.0   5.0)))

(defparameter *swg-aa-inserty*
  #2A((0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0)
      (0.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0)
      (0.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0 -10.0)
      (0.0 -10.0 -10.0  -5.0 -10.0  -5.0 -10.0 -10.0 -10.0 -10.0 -10.0)
      (0.0 -10.0 -10.0  -6.0  -8.0  -6.0  10.0   0.0  -1.0  -2.0  -3.0)
      (0.0   0.0 -10.0  -7.0  -9.0  -7.0   9.0   8.0  10.0   0.0  -1.0)
      (0.0  -1.0   6.0  -4.0  -5.0  -6.0   8.0   7.0   9.0  16.0   6.0)
      (0.0  -2.0   5.0  11.0   1.0   0.0   7.0   6.0   8.0  15.0  15.0)))

(defparameter *dna-query* (make-dna "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"))

(defun swg-aa-test (seqm seqn subst-fn gap-open gap-extend)
  (let ((vecm (bs::vector-of seqm))
        (vecn (bs::vector-of seqn)))
    (flet ((subn (x y) ; local fn to avoid boxing of returned floats
             (funcall subst-fn x y)))
      (let ((m (length vecm))
            (n (length vecn)))
        (bs::with-affine-gap-matrices
            (score insertx inserty backtrace) ((1+ m) (1+ n))
          (bs::define-affine-gap-dp
              (((row col) (prev-row prev-col) (max-row max-col))
               (cell-score max-score)
               (score insertx inserty backtrace (gap-open gap-extend))
               (subn (aref vecm prev-row) (aref vecn prev-col))
               ()) ; no cell exclusion form
            (values score insertx inserty))))))) ; return matrices for checking

(addtest (bio-sequence-alignment-tests) smith-waterman-gotoh-aa/1
  (let ((seqm (make-aa "PAWHEAE"))
        (seqn (make-aa "HEAGAWGHEE"))
        (subst-fn #'blosum-50-subst)
        (gap-open -10.0)
        (gap-extend -1.0))
    (multiple-value-bind (score insertx inserty)
        (swg-aa-test seqm seqn subst-fn gap-open gap-extend)
      (ensure (equalp *swg-aa-score* score))
      (ensure (equalp *swg-aa-insertx* insertx))
      (ensure (equalp *swg-aa-inserty* inserty)))
    (multiple-value-bind (align-score alignment)
        (align-local seqm seqn subst-fn :gap-open gap-open
                     :gap-extend gap-extend :alignment t)
      (ensure (= 26.0 align-score))
      (let ((intervalm (first (intervals-of alignment)))
            (intervaln (second (intervals-of alignment))))
        ;; Check coordinates of alignment
        (ensure (= 1 (lower-of intervalm)))
        (ensure (= 5 (upper-of intervalm))) ; interval like subseq
        (ensure (= 4 (lower-of intervaln)))
        (ensure (= 9 (upper-of intervaln))) ; interval like subseq
        ;; Check sequence and gapping
        (ensure (= 1 (num-gaps-of (aligned-of intervalm))))
        (ensure (= 0 (num-gaps-of (aligned-of intervaln))))
        (ensure (string= "AW-HE"
                         (coerce-sequence (aligned-of intervalm) 'string)))
        (ensure (string= "AWGHE"
                         (coerce-sequence (aligned-of intervaln) 'string)))))))

(addtest (bio-sequence-alignment-tests) smith-waterman-gotoh-dna/1
   (let ((subst-fn #'iupac-dna-subst)
         (gap-open -5.0)
         (gap-extend -1.0)
         ;; 1000 Fasta sequences
         (seqs (with-ascii-li-stream (stream (merge-pathnames
                                              "data/alignment_test_db.fasta"))
                 (loop
                    with gen = (make-seq-input stream :fasta)
                    as seq = (next gen)
                    while (has-more-p gen)
                    collect seq)))
         ;; Scores from aligning with EMBOSS 5.0 water
         (scores (with-open-file
                     (s (merge-pathnames "data/alignment-test-scores.sexp"))
                   (read s))))
     (mapc (lambda (seq score)
             (multiple-value-bind (s a)
                 (align-local *dna-query* seq subst-fn
                              :gap-open gap-open
                              :gap-extend gap-extend)
               (ensure (= score s))))
           seqs scores)))
