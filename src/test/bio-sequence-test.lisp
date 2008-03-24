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

;;; Test alphabet finding
(test find-alphabet/standard
  (is (eql *dna* (find-alphabet :dna)))
  (is (eql *rna* (find-alphabet :rna)))
  (is (eql *iupac-dna* (find-alphabet :iupac-dna)))
  (is (eql *iupac-rna* (find-alphabet :iupac-rna)))
  (signals error
    (find-alphabet :foo)))

;;; Test alphabet
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

(test standard-alphabets
  "Test presence of expected tokens in standard alphabets."
  (mapc #'(lambda (alphabet tokens)
            (is (eq (length tokens) (size-of alphabet)))
            (dolist (token tokens)
              (is-true (memberp alphabet token))))
        (mapcar #'find-alphabet '(:dna :rna :iupac-dna :iupac-rna))
        (mapcar #'(lambda (str)
                    (loop
                       for c across str collect c))
                '("acgt" "acgu" "acgtrykmswbdhvn" "acgurykmswbdhvn"))))

;;; Test encoding/decoding
(test encode/decode-dna-2bit
  (let ((tokens "tagc"))
    (loop for tok across tokens
         do (is (char= tok (bio-sequence::decode-dna-2bit
                            (bio-sequence::encode-dna-2bit tok)))))))

(test encode/decode-rna-2bit
  (let ((tokens "uagc"))
    (loop for tok across tokens
       do (is (char= tok (bio-sequence::decode-rna-2bit
                          (bio-sequence::encode-rna-2bit tok)))))))

(test encode/decode-dna-4bit
  (let ((tokens "tagcrykmswbdhvn"))
    (loop for tok across tokens
       do (is (char= tok (bio-sequence::decode-dna-4bit
                          (bio-sequence::encode-dna-4bit tok)))))))

(test encode/decode-rna-4bit
  (let ((tokens "uagcrykmswbdhvn"))
    (loop for tok across tokens
       do (is (char= tok (bio-sequence::decode-rna-4bit
                          (bio-sequence::encode-rna-4bit tok)))))))

(test encode-token
  (let* ((residues "acgt")
         (seq (make-instance 'simple-dna-sequence
                             :token-seq residues)))
    (loop for i from 0 below (length residues)
       do (progn
            (bio-sequence::encode-token (char residues i) seq i
                                        (bio-sequence::encoded-tokens 2))
            (is (char= (char residues i)
                       (char (to-string seq) i)))))))

(test dencode-token
  (let* ((residues "acgt")
         (seq (make-instance 'simple-dna-sequence
                             :token-seq residues)))
    (loop for i from 0 below (length residues)
       do (progn
            (is (char= (bio-sequence::decode-token
                        seq i
                        (bio-sequence::encoded-tokens 2))
                       (char residues i)))))))

;;; Test constructors
(test make-dna/rna
  (let ((seqs (list (make-seq :token-seq "tagc")
                    (make-seq :ambiguity :iupac :token-seq "nnnn")
                    (make-seq :alphabet :rna :token-seq "uagc")
                    (make-seq :alphabet :rna :ambiguity :iupac
                              :token-seq"nnnn")))
        (class-names (list 'simple-dna-sequence
                           'iupac-dna-sequence
                           'simple-rna-sequence
                           'iupac-rna-sequence)))
    (mapcar #'(lambda (seq class-name)
                (is (eql class-name (class-name (class-of seq)))))
            seqs class-names))
  (signals error
    (make-seq))
  (signals error
    (make-seq :token-seq ""))
  (signals error
    (make-seq :token-seq '(#\t #\a #\g #\c)))
  (signals error
    (make-seq :alphabet :dna :length -1))
  (signals error
    (make-seq :token-seq "tagc" :length 99)))

(test make-quality-dna
  (let ((seqs (list (make-quality-seq
                     :token-seq "agaatattctgaccccagttactttcaaga"
                     :quality "<<<<<<<<<<<<<<<<<<<<<735513;3<"
                     :metric :phred)
                    (make-quality-seq
                     :ambiguity :iupac
                     :token-seq "ntgccaaaaaatagaaaagtcancgatatt"
                     :quality "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
                     :metric :phred)))
        (class-names (list 'simple-dna-quality-sequence
                           'iupac-dna-quality-sequence)))
    (mapcar #'(lambda (seq class-name)
                (is (eql class-name (class-name (class-of seq)))))
            seqs class-names))
  (signals error
   (make-quality-seq
    :token-seq "agaatattctgaccccagttactttcaaga"
    :quality "<<<<<<<"
    :metric :phred))
  (signals error
    (make-quality-seq
     :token-seq "agaatattctgaccccagttactttcaaga"
     :quality "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
     :metric :invalid-metric)))

;; Test bio-sequence methods
(test simplep/string
  (is-true (simplep "acgt" :dna))
  (is-true (simplep "acgu" :rna)))

(test length-of/bio-sequence
  (let ((len 10))
    (is (= len (length-of (make-instance 'simple-dna-sequence
                                         :token-seq "aaaaaaaaaa"))))
    (is (= len (length-of (make-instance 'simple-dna-sequence
                                         :token-seq "aaaaaaaaaa")))))
  (signals error
    (setf (length-of (make-instance 'simple-dna-sequence
                                    :token-seq "aaaaaaaaaa")) 99)))

(test residue-of/dna-sequence
  (let ((residues "tttt")
        (seq (make-instance 'simple-dna-sequence
                            :token-seq "aaaa")))
    (dotimes (n (length residues))
      (setf (residue-of seq n) (aref residues n))
      (is (char= (residue-of seq n) (aref residues n))))))

(test residue-of/rna-sequence
  (let ((residues "uuuu")
        (rna-seq (make-instance 'simple-rna-sequence
                                :token-seq "aaaa")))
    (dotimes (n (length residues))
      (setf (residue-of rna-seq n) (aref residues n))
      (is (char= (residue-of rna-seq n) (aref residues n))))))

(test residue-of/iupac-dna-sequence
  (let ((residues "tntt")
        (seq (make-instance 'iupac-dna-sequence
                            :token-seq "aana")))
    (dotimes (n (length residues))
      (setf (residue-of seq n) (aref residues n))
      (is (char= (residue-of seq n) (aref residues n))))))

(test residue-of/iupac-rna-sequence
  (let ((residues "unuu")
        (rna-seq (make-instance 'iupac-rna-sequence
                                :token-seq "aana")))
    (dotimes (n (length residues))
      (setf (residue-of rna-seq n) (aref residues n))
      (is (char= (residue-of rna-seq n) (aref residues n))))))

(test residue-of/virtual-sequence
  (signals error
    (residue-of (make-instance 'simple-dna-sequence
                               :length 10) 0))
  (signals error
    (setf (residue-of (make-instance 'simple-dna-sequence
                                     :length 10) 0) #\a)))

;; (test copy-sequence/bio-sequence
;;   (let ((seqs (list (make-instance 'simple-dna-sequence
;;                                    :token-seq "tagc")
;;                     (make-instance 'iupac-dna-sequence
;;                                    :token-seq "nnnn")
;;                     (make-instance 'simple-rna-sequence
;;                                    :token-seq "uagc")
;;                     (make-instance 'iupac-rna-sequence
;;                                    :token-seq "nnnn"))))
;;     (mapcar #'(lambda (seq)
;;                 (let ((copy (copy-sequence seq)))
;;                   (is (= (length-of copy) (length-of seq)))
;;                   (dotimes (n 4)
;;                     (is (char= (residue-of copy n)
;;                                (residue-of seq n))))))
;;             seqs)))

(test to-string/dna-sequence
  (let* ((residues "acgt")
         (seq (make-instance 'simple-dna-sequence
                             :token-seq residues)))
    ;; no args
    (is (string= residues (to-string seq)))
    ;; optional arg start
     (dotimes (n 4)
       (is (string= (subseq residues n) (to-string seq n))))
     (dotimes (n 4)
       (is (string= (subseq residues 0 n) (to-string seq 0 n))))))

(test to-string/rna-sequence
  (let* ((residues "acgu")
         (seq (make-instance 'simple-rna-sequence
                             :token-seq residues)))
    ;; no args
    (is (string= residues (to-string seq)))
    ;; optional arg start
     (dotimes (n 4)
       (is (string= (subseq residues n) (to-string seq n))))
     (dotimes (n 4)
       (is (string= (subseq residues 0 n) (to-string seq 0 n))))))

(test subsequence/bio-sequence
  (let* ((residues "aaggccttaaggcctt")
         (seq (make-seq :token-seq residues)))
    (is (string= (subseq residues 0 5) ; aaggc
                 (to-string (subsequence seq 0 5))))
    (is (string= residues
                 (to-string (subsequence seq 0))))))

(test reverse-sequence/bio-sequence
  (let* ((residues "aaccggtt")
         (seq (make-seq :token-seq residues)))
    (is (string= (to-string (reverse-sequence seq))
                 (reverse residues)))))

(test nreverse-sequence/bio-sequence
  (let* ((residues "aaccggtt")
         (seq (make-seq :token-seq residues)))
    (is (string= (to-string (nreverse-sequence seq))
                 (reverse residues)))))

(test complement-sequence/simple-dna-sequence
  (let* ((residues "aaccggtt")
         (seq (make-seq :token-seq residues)))
    (is (string= (to-string (complement-sequence seq))
                 "ttggccaa"))))

(test complement-sequence/iupac-dna-sequence
  (let* ((residues "acgtrykmswbdhvn")
         (seq (make-seq :token-seq residues :ambiguity :iupac)))
    (is (string= (to-string (complement-sequence seq))
                 "tgcayrmkswvhdbn"))))

(test reverse-complement/simple-dna-sequence
  (let* ((residues "acccggttt")
         (seq (make-seq :token-seq residues)))
    (is (string= (to-string (reverse-complement seq))
                 "aaaccgggt"))))

(test reverse-complement/iupac-dna-sequence
  (let* ((residues "acgtrykmswbdhvn")
         (seq (make-seq :token-seq residues :ambiguity :iupac)))
    (is (string= (to-string (reverse-complement seq))
                 "nbdhvwskmryacgt"))))

(test residue-frequencies/bio-sequence
  (let ((seq (make-seq :token-seq "aacccgggt")))
    (mapc #'(lambda (token freq)
              (is (eq (cdr (assoc token (residue-frequencies seq)
                                  :test #'char= )) freq)))
          '(#\a #\c #\g #\t)
          '(2 3 3 1)))
  (signals error
    (residue-frequencies (make-seq :alphabet :dna :length 10))))
