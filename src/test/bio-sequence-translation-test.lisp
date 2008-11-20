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

(deftestsuite bio-sequence-translation-tests (cl-genomic-tests)
  ())

(defparameter *test-codons*
  '((#\a #\a #\a) (#\a #\a #\c) (#\a #\a #\g) (#\a #\a #\t) (#\a #\a #\r)
    (#\a #\a #\y) (#\a #\c #\a) (#\a #\c #\c) (#\a #\c #\g) (#\a #\c #\t)
    (#\a #\c #\r) (#\a #\c #\y) (#\a #\c #\k) (#\a #\c #\m) (#\a #\c #\s)
    (#\a #\c #\w) (#\a #\c #\b) (#\a #\c #\d) (#\a #\c #\h) (#\a #\c #\v)
    (#\a #\c #\n) (#\a #\g #\a) (#\a #\g #\c) (#\a #\g #\g) (#\a #\g #\t)
    (#\a #\g #\r) (#\a #\g #\y) (#\a #\t #\a) (#\a #\t #\c) (#\a #\t #\g)
    (#\a #\t #\t) (#\a #\t #\y) (#\a #\t #\m) (#\a #\t #\w) (#\a #\t #\h)
    (#\c #\a #\a) (#\c #\a #\c) (#\c #\a #\g) (#\c #\a #\t) (#\c #\a #\r)
    (#\c #\a #\y) (#\c #\c #\a) (#\c #\c #\c) (#\c #\c #\g) (#\c #\c #\t)
    (#\c #\c #\r) (#\c #\c #\y) (#\c #\c #\k) (#\c #\c #\m) (#\c #\c #\s)
    (#\c #\c #\w) (#\c #\c #\b) (#\c #\c #\d) (#\c #\c #\h) (#\c #\c #\v)
    (#\c #\c #\n) (#\c #\g #\a) (#\c #\g #\c) (#\c #\g #\g) (#\c #\g #\t)
    (#\c #\g #\r) (#\c #\g #\y) (#\c #\g #\k) (#\c #\g #\m) (#\c #\g #\s)
    (#\c #\g #\w) (#\c #\g #\b) (#\c #\g #\d) (#\c #\g #\h) (#\c #\g #\v)
    (#\c #\g #\n) (#\c #\t #\a) (#\c #\t #\c) (#\c #\t #\g) (#\c #\t #\t)
    (#\c #\t #\r) (#\c #\t #\y) (#\c #\t #\k) (#\c #\t #\m) (#\c #\t #\s)
    (#\c #\t #\w) (#\c #\t #\b) (#\c #\t #\d) (#\c #\t #\h) (#\c #\t #\v)
    (#\c #\t #\n) (#\g #\a #\a) (#\g #\a #\c) (#\g #\a #\g) (#\g #\a #\t)
    (#\g #\a #\r) (#\g #\a #\y) (#\g #\c #\a) (#\g #\c #\c) (#\g #\c #\g)
    (#\g #\c #\t) (#\g #\c #\r) (#\g #\c #\y) (#\g #\c #\k) (#\g #\c #\m)
    (#\g #\c #\s) (#\g #\c #\w) (#\g #\c #\b) (#\g #\c #\d) (#\g #\c #\h)
    (#\g #\c #\v) (#\g #\c #\n) (#\g #\g #\a) (#\g #\g #\c) (#\g #\g #\g)
    (#\g #\g #\t) (#\g #\g #\r) (#\g #\g #\y) (#\g #\g #\k) (#\g #\g #\m)
    (#\g #\g #\s) (#\g #\g #\w) (#\g #\g #\b) (#\g #\g #\d) (#\g #\g #\h)
    (#\g #\g #\v) (#\g #\g #\n) (#\g #\t #\a) (#\g #\t #\c) (#\g #\t #\g)
    (#\g #\t #\t) (#\g #\t #\r) (#\g #\t #\y) (#\g #\t #\k) (#\g #\t #\m)
    (#\g #\t #\s) (#\g #\t #\w) (#\g #\t #\b) (#\g #\t #\d) (#\g #\t #\h)
    (#\g #\t #\v) (#\g #\t #\n) (#\t #\a #\a) (#\t #\a #\c) (#\t #\a #\g)
    (#\t #\a #\t) (#\t #\a #\r) (#\t #\a #\y) (#\t #\c #\a) (#\t #\c #\c)
    (#\t #\c #\g) (#\t #\c #\t) (#\t #\c #\r) (#\t #\c #\y) (#\t #\c #\k)
    (#\t #\c #\m) (#\t #\c #\s) (#\t #\c #\w) (#\t #\c #\b) (#\t #\c #\d)
    (#\t #\c #\h) (#\t #\c #\v) (#\t #\c #\n) (#\t #\g #\a) (#\t #\g #\c)
    (#\t #\g #\g) (#\t #\g #\t) (#\t #\g #\y) (#\t #\t #\a) (#\t #\t #\c)
    (#\t #\t #\g) (#\t #\t #\t) (#\t #\t #\r) (#\t #\t #\y) (#\t #\r #\a)
    (#\r #\a #\c) (#\r #\a #\t) (#\r #\a #\y) (#\y #\t #\a) (#\y #\t #\g)
    (#\y #\t #\r) (#\m #\g #\a) (#\m #\g #\g) (#\m #\g #\r) (#\s #\a #\a)
    (#\s #\a #\g) (#\s #\a #\r))
  "Codons that do not encode for X.")

(defparameter *test-amino-acids*
  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVV*Y*Y*YSSSSSSSSSSSSSSS*CWCCLFLFLF*ZZZLLLRRRBBB"
  "Corresponding translations of codons that do not encode for X.")

(addtest (bio-sequence-translation-tests) translate-codon-standard
   (ensure (every #'char= *test-amino-acids*
                  (loop
                     with code = (find-genetic-code :standard)
                     for codon in *test-codons*
                     collect (bs::decode-aa-7bit
                              (translate-codon
                               (mapcar #'bs::encode-dna-4bit codon) code
                               :initiator nil)))))
   (ensure (every (lambda (char)
                    (or (char= #\X char)
                        (find char *test-amino-acids* :test #'char=)))
                  (loop
                     with code = (find-genetic-code :standard)
                     for codon in (bs::permutations-of-n
                                   '(#\a #\c #\g #\t #\r #\y #\k #\m
                                     #\s #\w #\b #\d #\h #\v #\n) 3)
                     collect (bs::decode-aa-7bit
                              (translate-codon
                               (mapcar #'bs::encode-dna-4bit codon) code
                               :initiator nil))))))
