;;;
;;; Copyright (C) 2008, Keith James. All rights reserved.
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

(in-package :cl-bio-test)

(fiveam:in-suite cl-bio-system:testsuite)

;;; Complementing nucleic acids
(test complement-dna
  (let ((residues "tagcrykmswbdhvn")
        (complemenent "atcgyrmkswvhdbn"))
    (loop
       for res across residues
       for cmp across complemenent
       do (progn
            (is (char= (bio-sequence::complement-dna res) cmp))
            (is (char= (bio-sequence::complement-dna (char-upcase res)) cmp))))))

(test complement-rna
  (let ((residues "uagcrykmswbdhvn")
        (complemenent "aucgyrmkswvhdbn"))
    (loop
       for res across residues
       for cmp across complemenent
       do (progn
            (is (char= (bio-sequence::complement-rna res) cmp))
            (is (char= (bio-sequence::complement-rna (char-upcase res)) cmp))))))

(test complement-dna-4bit
  (let ((residues     "tagcrykmswbdhvn")
        (complemenent "atcgyrmkswvhdbn"))
    (loop
       for res across residues
       for cmp across complemenent
       do (progn
            (is (char= (bio-sequence::decode-dna-4bit
                        (bio-sequence::complement-dna-4bit
                         (bio-sequence::encode-dna-4bit res)))
                       cmp))))))

;;; Encoding/decoding sequences
(test encode/decode-dna-4bit
  (let ((residues "tagcrykmswbdhvn"))
    (loop
       for res across residues
       do (is (char= res (bio-sequence::decode-dna-4bit
                          (bio-sequence::encode-dna-4bit res)))))))

(test encode/decode-rna-4bit
  (let ((residues "uagcrykmswbdhvn"))
    (loop
       for res across residues
       do (is (char= res (bio-sequence::decode-rna-4bit
                          (bio-sequence::encode-rna-4bit res)))))))

(test explode-dna-4bit
  (let ((residues "tagcrykmswbdhvn")
        (exploded (pairlis '(#\t #\a #\g #\c
                             #\r #\y #\k #\m
                             #\s #\w #\b #\d
                             #\h #\v #\n)
                           '((#\t) (#\a) (#\g) (#\c)
                             (#\a #\g) (#\c #\t) (#\g #\t) (#\a #\c)
                             (#\c #\g) (#\a #\t) (#\c #\g #\t) (#\a #\g #\t)
                             (#\a #\c #\t) (#\a #\c #\g) (#\a #\c #\g #\t)))))
    (loop
       for res across residues
       do (is (equal (assocdr res exploded :test #'char=)
                     (explode-ambiguity (find-alphabet :dna) res))))))

(test ambiguousp
  (let ((unambiguous "tagc")
        (ambiguous "rykmswbdhvn"))
    (loop
       for res across unambiguous
       do (is-false (ambiguousp (find-alphabet :dna) res)))
    (loop
       for res across ambiguous
       do (is-true (ambiguousp (find-alphabet :dna) res)))))
