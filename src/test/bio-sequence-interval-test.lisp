;;;
;;; Copyright (C) 2008 Keith. All rights reserved.
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

(deftestsuite bio-sequence-interval-tests (cl-genomic-tests)
  ((ss-seq (make-dna dna-residues))
   (ds-seq (make-dna dna-residues :num-strands 2))))

(addtest (bio-sequence-interval-tests) interval/1
  (let* ((lower 2)
         (upper 5)
         (interval (make-instance 'interval
                                  :reference ss-seq
                                  :lower lower
                                  :upper upper)))
    (ensure (= lower (lower-of interval)))
    (ensure (= upper (upper-of interval)))
    (ensure-same ss-seq (reference-of interval))
    (ensure-error 'invalid-argument-error
                  (setf (reference-of interval) (make-dna "a")))
    (ensure-error 'invalid-argument-error
                  (setf (lower-of interval) -1))
    (ensure-error 'invalid-argument-error
                  (setf (lower-of interval) 6))
    (ensure-error 'invalid-argument-error
                  (setf (upper-of interval) -1))
    (ensure-error 'invalid-argument-error
                  (setf (upper-of interval) 1))))

(addtest (bio-sequence-interval-tests) bio-sequence-interval/1
  (ensure-error 'invalid-argument-error
                (make-instance 'bio-sequence-interval
                               :lower 10 :upper 1)))

(addtest (bio-sequence-interval-tests) bio-sequence-interval/2
  (ensure (= 5 (length-of (make-instance 'bio-sequence-interval
                                         :reference ss-seq
                                         :lower 0
                                         :upper 5)))))

(addtest (bio-sequence-interval-tests) bio-sequence-interval/3
  (ensure (string= "-----"
                   (coerce-sequence
                    (make-instance 'bio-sequence-interval
                                   :reference nil
                                   :lower 0
                                   :upper 5) 'string))))

(addtest (bio-sequence-interval-tests) bio-sequence-interval/4
  (ensure (string= "tagcr"
                   (coerce-sequence
                    (make-instance 'bio-sequence-interval
                                   :reference ss-seq
                                   :lower 0
                                   :upper 5) 'string))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/1
  (ensure-error 'invalid-argument-error
                (make-instance 'na-sequence-interval
                               :lower 10 :upper 1))
  (ensure-error 'invalid-argument-error
                (make-instance 'na-sequence-interval
                               :reference (make-dna "acgt" :num-strands 1)
                               :lower 0
                               :upper 5
                               :strand *reverse-strand*))
  (ensure-error 'invalid-argument-error
                (make-instance 'na-sequence-interval
                               :lower 0
                               :upper 5
                               :strand *forward-strand*
                               :reference (make-dna "acgt" :num-strands 1)
                               :num-strands 2)
                :report "Made a double-stranded interval with a
single-stranded reference."))

(addtest (bio-sequence-interval-tests) na-sequence-interval/2
  (let* ((lower 2)
         (upper 5)
         (interval1 (make-instance 'na-sequence-interval
                                   :reference ss-seq
                                   :lower lower
                                   :upper upper))
         (interval2 (make-instance 'na-sequence-interval
                                   :reference ds-seq
                                   :lower lower
                                   :upper upper
                                   :strand *forward-strand*))
         (interval3 (make-instance 'na-sequence-interval
                                   :reference ds-seq
                                   :lower lower
                                   :upper upper
                                   :strand *forward-strand*
                                   :num-strands 2))
         (interval4 (make-instance 'na-sequence-interval
                                   :reference ds-seq
                                   :lower lower
                                   :upper upper
                                   :strand *reverse-strand*)))
    (ensure (= lower (lower-of interval1)))
    (ensure (= upper (upper-of interval1)))
    (ensure-same ss-seq (reference-of interval1))
    (ensure-error 'invalid-argument-error
                  (setf (reference-of interval3) ss-seq))
    (ensure-error 'invalid-argument-error
                  (setf (reference-of interval4) ss-seq))
    (ensure-error 'invalid-argument-error
                  (setf (num-strands-of interval1) 2))
    (ensure-error 'invalid-argument-error
                  (setf (strand-of interval1) *reverse-strand*))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/3
  (ensure (= 5 (length-of (make-instance 'na-sequence-interval
                                         :reference ss-seq
                                         :lower 0
                                         :upper 5)))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/4
  (ensure (string= "-----"
                   (coerce-sequence
                    (make-instance 'na-sequence-interval
                                   :reference nil
                                   :lower 0
                                   :upper 5) 'string))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/5
  (ensure (string= "tagcr"
                   (coerce-sequence
                    (make-instance 'na-sequence-interval
                                   :reference ss-seq
                                   :lower 0
                                   :upper 5
                                   :strand *forward-strand*) 'string))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/6
  (ensure (string= "ygcta"
                   (coerce-sequence
                    (make-instance 'na-sequence-interval
                                   :reference ds-seq
                                   :lower 0
                                   :upper 5
                                   :strand *reverse-strand*) 'string))))

(let* ((reference1 (make-dna "aagggct"))
       (reference2 (make-dna "aagggct" :num-strands 2)))

  (addtest (bio-sequence-interval-tests) na-sequence-interval/7
    (let ((interval1 (make-instance 'na-sequence-interval
                                     :reference reference1
                                     :lower 1
                                     :upper 6)))
      ;; interval *unknown-strand*
      (ensure (string= "aggg"
                       (coerce-sequence (subsequence interval1 0 4) 'string)))
      (ensure (string= "ggg"
                       (coerce-sequence (subsequence interval1 1 4) 'string)))
      (ensure (string= "gggc"
                       (coerce-sequence (subsequence interval1 1) 'string)))))

  (addtest (bio-sequence-interval-tests) na-sequence-interval/8
    (let ((interval2 (make-instance 'na-sequence-interval
                                    :reference reference1
                                    :lower 1
                                    :upper 6
                                    :strand *forward-strand*)))
      ;; interval *forward-strand*
      (ensure (string= "aggg"
                       (coerce-sequence (subsequence interval2 0 4) 'string)))
      (ensure (string= "ggg"
                       (coerce-sequence (subsequence interval2 1 4) 'string)))
      (ensure (string= "gggc"
                       (coerce-sequence (subsequence interval2 1) 'string)))))

  (addtest (bio-sequence-interval-tests) na-sequence-interval/9
    (let ((interval3 (make-instance 'na-sequence-interval
                                    :reference reference2
                                    :lower 1
                                    :upper 6
                                    :strand *reverse-strand*)))
      ;; interval *reverse-strand*
      (ensure (string= "ccct"
                       (coerce-sequence (subsequence interval3 0 4) 'string)))
      (ensure (string= "ccc"
                       (coerce-sequence (subsequence interval3 1 4) 'string)))
      (ensure (string= "gccc"
                       (coerce-sequence (subsequence interval3 1) 'string))))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/10
  (let ((reference-ss (make-dna "aaacgtttgc"))
        (reference-ds (make-dna "aaacgtttgc" :num-strands 2))
        (lower 0)
        (upper 4)
        (strand *forward-strand*))
    (let* ((interval1 (make-instance 'na-sequence-interval
                                     :lower lower
                                     :upper upper
                                     :strand strand
                                     :reference reference-ds))
           (interval2 (invert-complement interval1))
           (interval3 (invert-complement interval2))
           (interval4 (make-instance 'na-sequence-interval
                                     :lower lower
                                     :upper upper
                                     :strand strand
                                     :reference reference-ss)))
      ;; interval2 is the reverse-complement of interval1
      (ensure (eql *reverse-strand* (strand-of interval2)))
      (ensure (= 6 (lower-of interval2)))
      (ensure (= 10 (upper-of interval2)))
      (ensure (string= "gcaa" (coerce-sequence interval2 'string)))
      ;; interval3 is the reverse-complement of interval2, so should
      ;; be identical to interval1
      (ensure (eql *forward-strand* (strand-of interval3)))
      (ensure (= 0 (lower-of interval3)))
      (ensure (= 4 (upper-of interval3)))
      (ensure (string= "aaac" (coerce-sequence interval3 'string)))
      (ensure-error 'invalid-argument-error
                    (invert-complement interval4)
                    :report "Successfully called invert-complement on
                    an interval with a single-stranded reference."))))

(addtest (bio-sequence-interval-tests) na-sequence-interval/11
  (let ((reference-ss (make-dna "aaacgtttgc"))
        (reference-ds (make-dna "aaacgtttgc" :num-strands 2))
        (lower 0)
        (upper 4)
        (strand *forward-strand*))
    (let* ((interval1 (make-instance 'na-sequence-interval
                                     :lower lower
                                     :upper upper
                                     :strand strand
                                     :reference reference-ds))
           (interval2 (make-instance 'na-sequence-interval
                                     :lower lower
                                     :upper upper
                                     :strand strand
                                     :reference reference-ss)))
      (ninvert-complement interval1)
      (ensure (eql *reverse-strand* (strand-of interval1)))
      (ensure (= 6 (lower-of interval1)))
      (ensure (= 10 (upper-of interval1)))
      (ensure (string= "gcaa" (coerce-sequence interval1 'string)))
      (ninvert-complement interval1)
      (ensure (eql *forward-strand* (strand-of interval1)))
      (ensure (= 0 (lower-of interval1)))
      (ensure (= 4 (upper-of interval1)))
      (ensure (string= "aaac" (coerce-sequence interval1 'string)))
      (ensure-error 'invalid-argument-error
                    (invert-complement interval2)
                    :report "Successfully called ninvert-complement on
                    an interval with a single-stranded reference."))))
