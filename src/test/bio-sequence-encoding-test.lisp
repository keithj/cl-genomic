;;;
;;; Copyright (C) 2008-2009 Keith James. All rights reserved.
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

(deftestsuite bio-sequence-encoding-tests (cl-genomic-tests)
  ((dna-complement"atcgyrmkswvhdbn-")
   (rna-complement"aucgyrmkswvhdbn-")
   (dna-unambiguous "tagc")
   (dna-ambiguous "rykmswbdhvn")
   (dna-enum (pairlis '(#\t #\a #\g #\c
                        #\r #\y #\k #\m
                        #\s #\w #\b #\d
                        #\h #\v #\n)
                      '((#\t) (#\a) (#\g) (#\c)
                        (#\a #\g) (#\c #\t) (#\g #\t) (#\a #\c)
                        (#\c #\g) (#\a #\t) (#\c #\g #\t) (#\a #\g #\t)
                        (#\a #\c #\t) (#\a #\c #\g) (#\a #\c #\g #\t))))))

;;; Complementing nucleic acids
(addtest (bio-sequence-encoding-tests) complement-dna/1
  (loop
     for res across *dna-residues*
     for cmp across dna-complement
     do (ensure (char= (bs::complement-dna res) cmp))))

(addtest (bio-sequence-encoding-tests) complement-dna/2
  (loop
     for res across *dna-residues*
     for cmp across dna-complement
     do (ensure (char= (bs::complement-dna (char-upcase res)) cmp))))

(addtest (bio-sequence-encoding-tests) complement-rna/1
  (loop
     for res across *rna-residues*
     for cmp across rna-complement
     do (ensure (char= (bs::complement-rna res) cmp))))

(addtest (bio-sequence-encoding-tests) complement-rna/2
  (loop
     for res across *rna-residues*
     for cmp across rna-complement
     do (ensure (char= (bs::complement-rna (char-upcase res)) cmp))))

(addtest (bio-sequence-encoding-tests) complement-dna-4bit/1
  (loop
     for res across *dna-residues*
     for cmp across dna-complement
     do (progn
          (ensure (char= (decode-dna-4bit
                          (bs::complement-dna-4bit
                           (encode-dna-4bit res)))
                         cmp)))))

;;; Encoding/decoding sequences
(addtest (bio-sequence-encoding-tests) encode/decode-dna-4bit/1
  (loop
     for res across *dna-residues*
     do (ensure (char= res (decode-dna-4bit (encode-dna-4bit res))))))

(addtest (bio-sequence-encoding-tests) encode/decode-rna-4bit/1
  (loop
     for res across *rna-residues*
     do (ensure (char= res (decode-rna-4bit (encode-rna-4bit res))))))

(addtest (bio-sequence-encoding-tests) enum-ambiguity/1
  (loop
     for res across *dna-residues*
     do (ensure (equal (assocdr res dna-enum :test #'char=)
                       (enum-ambiguity res (find-alphabet :dna))))))
