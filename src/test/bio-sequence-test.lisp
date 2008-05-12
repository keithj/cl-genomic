;;;
;;; Copyright (C) 2007-2008, Keith James. All rights reserved.
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

(in-package :cl-bio-system)

(fiveam:def-suite testsuite
    :description "The test suite.")


(in-package :cl-bio-test)

(fiveam:in-suite cl-bio-system:testsuite)

;;; Alphabets
(test find-alphabet/standard
  (is (eql *dna* (find-alphabet :dna)))
  (is (eql *rna* (find-alphabet :rna)))
  (signals error
    (find-alphabet :foo)))

(test standard-alphabets
  "Test presence of expected tokens in standard alphabets."
  (mapc #'(lambda (alphabet tokens)
            (is (eq (length tokens) (size-of alphabet)))
            (dolist (token tokens)
              (is-true (memberp alphabet token))))
        (mapcar #'find-alphabet '(:dna :rna))
        (mapcar #'(lambda (str)
                    (loop
                       for c across str collect c))
                '("acgtrykmswbdhvn" "acgurykmswbdhvn"))))

(test name-of/alphabet
  "Test alphabet."
  (let ((name "test name"))
    (is (string= name
                 (name-of (make-instance 'alphabet
                                         :name name))))))

(test tokens-of/alphabet
  (let ((tokens "abcd"))
    (is (string= tokens
                 (tokens-of (make-instance 'alphabet
                                           :tokens tokens))))))

;;; Encoding/decoding sequences
(test encode/decode-dna-4bit
  (let ((residues "tagcrykmswbdhvn"))
    (loop for res across residues
       do (is (char= res (bio-sequence::decode-dna-4bit
                          (bio-sequence::encode-dna-4bit res)))))))

(test encode/decode-rna-4bit
  (let ((residues "uagcrykmswbdhvn"))
    (loop for res across residues
       do (is (char= res (bio-sequence::decode-rna-4bit
                          (bio-sequence::encode-rna-4bit res)))))))

;;; Sequence constructors
(test make-dna/rna
  (let ((seqs (list (make-instance 'dna-sequence
                                   :residues "tagc") ; unambiguous
                    (make-instance 'dna-sequence
                                   :residues "nnnn") ; ambiguous
                    (make-instance 'rna-sequence
                                   :residues "uagc") ; unambiguous
                    (make-instance 'rna-sequence
                                   :residues "nnnn"))) ; ambiguous
        (class-names (list 'dna-sequence
                           'dna-sequence
                           'rna-sequence
                           'rna-sequence))
        (alphabets (list *dna*
                         *dna*
                         *rna*
                         *rna*)))
    (mapcar #'(lambda (seq class-name alphabet)
                (is (eql class-name (class-name (class-of seq))))
                (is (eql alphabet (alphabet-of seq))))
            seqs class-names alphabets))
  (signals error
    (make-instance 'dna-sequence))
  (signals error
    (make-instance 'dna-sequence :residues "u"))
  (signals error
    (make-instance 'rna-sequence :residues "t"))
  (signals error
    (make-instance 'dna-sequence :residues ""))
  (signals error
    (make-instance 'dna-sequence :residues '(#\t #\a #\g #\c)))
  (signals error
    (make-instance 'dna-sequence :length -1))
  (signals error
    (make-instance 'dna-sequence :residues "tagc" :length 99)))

(test make-quality-dna
  (let ((seqs (list (make-instance 'dna-quality-sequence ; unambiguous
                     :residues "agaatattctgaccccagttactttcaaga"
                     :quality "<<<<<<<<<<<<<<<<<<<<<735513;3<"
                     :metric :phred)
                    (make-instance 'dna-quality-sequence ; ambiguous
                     :residues "ntgccaaaaaatagaaaagtcancgatatt"
                     :quality "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
                     :metric :phred)))
        (class-names (list 'dna-quality-sequence
                           'dna-quality-sequence)))
    (mapcar #'(lambda (seq class-name)
                (is (eql class-name (class-name (class-of seq)))))
            seqs class-names))
  (signals error
   (make-instance 'dna-quality-sequence
    :residues "agaatattctgaccccagttactttcaaga"
    :quality "<<<<<<<"
    :metric :phred))
  (signals error
    (make-instance 'dna-quality-sequence
     :residues "agaatattctgaccccagttactttcaaga"
     :quality "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
     :metric :invalid-metric)))

;;; Utility methods
(test simplep/string
  (is-true (simplep "acgt" (find-alphabet :dna)))
  (is-true (simplep "acgu" (find-alphabet :rna))))


;;; Sequence accessors
(test length-of/bio-sequence
  (let ((len 10))
    (is (= len (length-of (make-instance 'dna-sequence
                                         :residues "aaaaaaaaaa"))))
    (is (= len (length-of (make-instance 'dna-sequence
                                         :residues "aaaaaaaaaa")))))
  (signals error
    (setf (length-of (make-instance 'dna-sequence
                                    :residues "aaaaaaaaaa")) 99)))

(test residue-of/dna-sequence
  (let ((residues "tttt")
        (seq (make-instance 'dna-sequence
                            :residues "aaaa")))
    (dotimes (n (length residues))
      (setf (residue-of seq n) (aref residues n))
      (is (char= (residue-of seq n) (aref residues n))))))

(test residue-of/rna-sequence
  (let ((residues "uuuu")
        (rna-seq (make-instance 'rna-sequence
                                :residues "aaaa")))
    (dotimes (n (length residues))
      (setf (residue-of rna-seq n) (aref residues n))
      (is (char= (residue-of rna-seq n) (aref residues n))))))

(test residue-of/virtual-sequence
  (signals error
    (residue-of (make-instance 'dna-sequence
                               :length 10) 0))
  (signals error
    (setf (residue-of (make-instance 'dna-sequence
                                     :length 10) 0) #\a)))

;;; Sequence transformations
(test to-string/dna-sequence
  (let* ((residues "acgt")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    ;; no args
    (is (string= residues (to-string seq)))
    ;; optional arg start
     (dotimes (n 4)
       (is (string= (subseq residues n) (to-string seq :start n))))
     (dotimes (n 4)
       (is (string= (subseq residues 0 n) (to-string seq :start 0 :end n))))))

(test to-string/rna-sequence
  (let* ((residues "acgu")
         (seq (make-instance 'rna-sequence
                             :residues residues)))
    ;; no args
    (is (string= residues (to-string seq)))
    ;; optional arg start
     (dotimes (n 4)
       (is (string= (subseq residues n) (to-string seq :start n))))
     (dotimes (n 4)
       (is (string= (subseq residues 0 n) (to-string seq :start 0 :end n))))))

(test reverse-sequence/dna-sequence
  (let* ((residues "aaccggtt")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    (is (string= (to-string (reverse-sequence seq))
                 (reverse residues)))))

(test reverse-sequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (seq (make-instance 'dna-quality-sequence 
                             :residues residues
                             :quality quality
                             :metric :phred))
         (rseq (reverse-sequence seq)))
    (is (string= (to-string rseq)
                 (reverse residues)))
    (is (string= (reverse quality)
                 (map-into (make-string (length residues))
                           #'encode-phred-quality
                           (quality-of rseq))))))

(test nreverse-sequence/dna-sequence
  (let* ((residues "aaccggtt")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    (is (string= (to-string (nreverse-sequence seq))
                 (reverse residues)))))

(test nreverse-sequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (rresidues (reverse residues))
         (rquality (reverse quality))
         (seq (make-instance 'dna-quality-sequence 
                             :residues residues
                             :quality quality
                             :metric :phred))
         (rseq (nreverse-sequence seq)))
    (is (string= rresidues (to-string rseq)))
    (is (equalp rquality
                (map-into (make-string (length residues))
                          #'encode-phred-quality
                          (quality-of rseq))))))

(test complement-sequence/dna-sequence
  (let* ((residues "acgtrykmswbdhvn")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    (is (string= (to-string (complement-sequence seq))
                 "tgcayrmkswvhdbn"))))

(test complement-sequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (seq (make-instance 'dna-quality-sequence 
                             :residues residues
                             :quality quality
                             :metric :phred)))
    (is (string= "tcttataagactggggtcaatgaaagttct"
                 (to-string (complement-sequence seq))))
    (is (string= quality
                 (map-into (make-string (length residues))
                          #'encode-phred-quality
                          (quality-of seq))))))

(test reverse-complement/dna-sequence
  (let* ((residues "acgtrykmswbdhvn")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    (is (string= (to-string (reverse-complement seq))
                 "nbdhvwskmryacgt"))))

(test reverse-complement/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (rquality (reverse quality))
         (seq (make-instance 'dna-quality-sequence 
                             :residues residues
                             :quality quality
                             :metric :phred)))
    (is (string= "tcttgaaagtaactggggtcagaatattct"
                 (to-string (reverse-complement seq))))
    (let ((rcquality (quality-of (reverse-complement seq))))
      (is (string= rquality
                   (map-into (make-string (length residues))
                             #'encode-phred-quality rcquality))))))

(test nreverse-complement/dna-sequence
  (let* ((residues "acgtrykmswbdhvn")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    (is (string= (to-string (nreverse-complement seq))
                 "nbdhvwskmryacgt"))))

(test nreverse-complement/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (rquality (reverse quality))
         (seq (make-instance 'dna-quality-sequence 
                             :residues residues
                             :quality quality
                             :metric :phred))
         (rseq (nreverse-complement seq)))
    (is (string= "tcttgaaagtaactggggtcagaatattct" (to-string rseq)))
    (let ((rcquality (quality-of rseq)))
      (is (string= rquality
                   (map-into (make-string (length residues))
                             #'encode-phred-quality rcquality))))))

(test subsequence/bio-sequence
  (let* ((residues "aaggccttaaggcctt")
         (seq (make-instance 'dna-sequence
                             :residues residues)))
    (is (string= (subseq residues 0 5) ; aaggc
                 (to-string (subsequence seq 0 5))))
    (is (string= residues
                 (to-string (subsequence seq 0))))))

(test subsequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (seq (make-instance 'dna-quality-sequence
                             :residues residues
                             :quality quality
                             :metric :phred)))
    (is (string= (subseq residues 0 5)
                 (to-string (subsequence seq 0 5))))
    (loop
       for q across (quality-of (subsequence seq 0 5))
       do (is (= 27 q)))))


(test residue-frequencies/bio-sequence
  (let ((seq (make-instance 'dna-sequence
                            :residues "aacccgggt")))
    (mapc #'(lambda (token freq)
              (is (eq (cdr (assoc token (residue-frequencies seq)
                                  :test #'char= )) freq)))
          '(#\a #\c #\g #\t)
          '(2 3 3 1)))
  (signals error
    (residue-frequencies (make-instance 'dna-sequence
                                        :length 10))))
