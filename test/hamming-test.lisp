;;;
;;; Copyright (c) 2010-2011 Keith James. All rights reserved.
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

(in-package :cl-genomic-test)

(deftestsuite hamming-tests (cl-genomic-tests)
  ())

(defun hamming-search-test (seq1 seq2 position distance &rest args)
  (multiple-value-bind (pos dist)
      (apply #'hamming-search seq1 seq2 args)
    (ensure (eql position pos))
    (ensure (eql distance dist))))

(addtest (hamming-tests) hamming-distance/1
  (ensure (= 0 (hamming-distance (make-dna "a") (make-dna "a"))))
  (ensure (= 1 (hamming-distance (make-dna "a") (make-dna "t"))))
  (ensure (= 0 (hamming-distance (make-dna "acgt")  (make-dna "acgt"))))

  ;; Hamming distance does not use base ambiguities
  (ensure (= 1 (hamming-distance (make-dna "acgt")  (make-dna "ncgt"))))
  (ensure (= 2 (hamming-distance (make-dna "acgt")  (make-dna "nngt"))))
  (ensure (= 3 (hamming-distance (make-dna "acgt")  (make-dna "nnnt"))))
  (ensure (= 4 (hamming-distance (make-dna "acgt")  (make-dna "nnnn"))))

  (ensure (= 0 (hamming-distance (make-dna "acgt")  (make-dna "acgtt"))))
  (ensure (= 0 (hamming-distance (make-dna "acgtt") (make-dna "acgt"))))
  (ensure (= 3 (hamming-distance (make-dna "aacgt") (make-dna "acgt"))))
  (ensure (= 3 (hamming-distance (make-dna "acgt")  (make-dna "aacgt")))))

(addtest (hamming-tests) hamming-distance/2
  (ensure (= 0 (hamming-distance (make-dna "aacgt") (make-dna "acgt")
                                 :start1 1)))
  (ensure (= 0 (hamming-distance (make-dna "acgt") (make-dna "nacgt")
                                 :start2 1)))
  (ensure (= 0 (hamming-distance (make-dna "acgtn") (make-dna "acgtt")
                                 :end1 3)))
  (ensure (= 0 (hamming-distance (make-dna "acgtt") (make-dna "acgtn")
                                 :end2 3)))
  (ensure (= 0 (hamming-distance (make-dna "aattgg") (make-dna "ggttcc")
                                 :start1 2 :end1 4 :start2 2 :end2 4))))

(addtest (hamming-tests) hamming-distance/3
  (ensure-error
    (hamming-distance (make-dna "") (make-dna ""))))

(addtest (hamming-tests) hamming-search/1
  (hamming-search-test (make-dna "acgt") (make-dna "acgt") 0 0)

  (hamming-search-test (make-dna "acgt") (make-dna "acgtnnnn") 0 0)
  (hamming-search-test (make-dna "acgt") (make-dna "nacgtnnn") 1 0)
  (hamming-search-test (make-dna "acgt") (make-dna "nnacgtnn") 2 0)
  (hamming-search-test (make-dna "acgt") (make-dna "nnnacgtn") 3 0)
  (hamming-search-test (make-dna "acgt") (make-dna "nnnnacgt") 4 0))

(addtest (hamming-tests) hamming-search/2
  ;; max-distance defaults to 1
  (hamming-search-test (make-dna "ayyt") (make-dna "acgtnnnn") nil nil)
  (hamming-search-test (make-dna "ayyt") (make-dna "nacgtnnn") nil nil)
  (hamming-search-test (make-dna "ayyt") (make-dna "nnacgtnn") nil nil)
  (hamming-search-test (make-dna "ayyt") (make-dna "nnnacgtn") nil nil)
  (hamming-search-test (make-dna "ayyt") (make-dna "nnnnacgt") nil nil))

(addtest (hamming-tests) hamming-search/3
  (hamming-search-test (make-dna "ayyt") (make-dna "acgtnnnn") 0 2
                       :max-distance 2)
  (hamming-search-test (make-dna "ayyt") (make-dna "nacgtnnn") 1 2
                       :max-distance 2)
  (hamming-search-test (make-dna "ayyt") (make-dna "nnacgtnn") 2 2
                       :max-distance 2)
  (hamming-search-test (make-dna "ayyt") (make-dna "nnnacgtn") 3 2
                       :max-distance 2)
  (hamming-search-test (make-dna "ayyt") (make-dna "nnnnacgt") 4 2
                       :max-distance 2))

(addtest (hamming-tests) hamming-search/4
  (hamming-search-test (make-dna "acgtnnnn") (make-dna "acgt") nil nil)
  (hamming-search-test (make-dna "acgtnnnn") (make-dna "acgt") 0 0
                       :end1 4)
  (hamming-search-test (make-dna "nnnnacgtnnnn") (make-dna "acgt") 0 0
                        :start1 4 :end1 8)
  (hamming-search-test (make-dna "nnnnacgtnnnn") (make-dna "yyyyacgt") 4 0
                       :start1 4 :end1 8 :start2 4)
  (hamming-search-test (make-dna "nnnnacgtnnnn") (make-dna "yyyyacgt") nil nil
                       :start1 4 :end1 8 :start2 5)
  (hamming-search-test (make-dna "nnnnacgtnnnn") (make-dna "yyyyacgtr") 5 4
                       :start1 4 :end1 8 :start2 5 :max-distance 4))
