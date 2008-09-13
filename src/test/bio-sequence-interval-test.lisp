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

(addtest (bio-sequence-interval-tests)
  make-instance/bio-sequence-interval
  (ensure-error 'invalid-argument-error
                (make-instance 'bio-sequence-interval
                               :lower 10 :upper 1)))

(addtest (bio-sequence-interval-tests)
  make-instance/nucleic-acid-sequence-interval
  (ensure-error 'invalid-argument-error
                (make-instance 'nucleic-acid-sequence-interval
                               :lower 10 :upper 1))
  (ensure-error 'invalid-argument-error
                (make-instance 'nucleic-acid-sequence-interval
                               :reference (make-dna "acgt" :num-strands 1)
                               :lower 0
                               :upper 5
                               :strand *reverse-strand*)))

(addtest (bio-sequence-interval-tests) bio-sequence-interval/accessors
  (let* ((lower 0)
         (upper 5)
         (i (make-instance 'bio-sequence-interval
                           :reference ss-seq
                           :lower lower
                           :upper upper)))
    (ensure (= lower (lower-of i)))
    (ensure (= upper (upper-of i)))
    (ensure-same ss-seq (reference-of i))))

(addtest (bio-sequence-interval-tests)
  nucleic-acid-sequence-interval/accessors
  (let* ((lower 0)
         (upper 5)
         (i (make-instance 'nucleic-acid-sequence-interval
                           :reference ss-seq
                           :lower lower
                           :upper upper)))
    (ensure (= lower (lower-of i)))
    (ensure (= upper (upper-of i)))
    (ensure-same ss-seq (reference-of i))))

(addtest (bio-sequence-interval-tests) length-of/bio-sequence-interval
  (ensure (= 5 (length-of (make-instance 'bio-sequence-interval
                                         :reference ss-seq
                                         :lower 0
                                         :upper 5)))))

(addtest (bio-sequence-interval-tests) length-of/nucleic-acid-sequence-interval
  (ensure (= 5 (length-of (make-instance 'nucleic-acid-sequence-interval
                                         :reference ss-seq
                                         :lower 0
                                         :upper 5)))))

(addtest (bio-sequence-interval-tests) to-string/bio-sequence-interval
  (ensure (string= ""
                   (to-string (make-instance 'bio-sequence-interval
                                             :reference nil
                                             :lower 0
                                             :upper 5))))
  (ensure (string= "tagcr"
                   (to-string (make-instance 'bio-sequence-interval
                                             :reference ss-seq
                                             :lower 0
                                             :upper 5)))))

(addtest (bio-sequence-interval-tests) to-string/nucleic-acid-sequence-interval
  (ensure (string= ""
                   (to-string (make-instance 'nucleic-acid-sequence-interval
                                             :reference nil
                                             :lower 0
                                             :upper 5))))
  (ensure (string= "tagcr"
                   (to-string
                    (make-instance 'nucleic-acid-sequence-interval
                                   :reference ss-seq
                                   :lower 0
                                   :upper 5
                                   :strand *forward-strand*))))
  (ensure (string= "atcgy"
                   (to-string
                    (make-instance 'nucleic-acid-sequence-interval
                                   :reference ds-seq
                                   :lower 0
                                   :upper 5
                                   :strand *reverse-strand*)))))
