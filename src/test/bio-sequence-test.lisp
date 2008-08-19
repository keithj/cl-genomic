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

(in-package :cl-genomic-test)

(deftestsuite bio-sequence-tests (cl-genomic-tests)
  ())

;;; Alphabets
(addtest (bio-sequence-tests) find-alphabet/standard
  (ensure-same *dna* (find-alphabet :dna))
  (ensure-same *rna* (find-alphabet :rna))
  (ensure-error
    (find-alphabet :foo)))

(addtest (bio-sequence-tests) standard-alphabets
  (mapc #'(lambda (alphabet tokens)
            (ensure (eq (length tokens) (size-of alphabet)))
            (ensure
             (loop
                for token in tokens
                always (memberp alphabet token))))
        (mapcar #'find-alphabet '(:dna :rna))
        (mapcar #'(lambda (str)
                    (loop
                       for c across str collect c))
                (list dna-residues rna-residues))))

(addtest (bio-sequence-tests) name-of/alphabet
  (let ((name "test name"))
    (ensure (string= name
                     (name-of (make-instance 'alphabet
                                             :name name))))))

(addtest (bio-sequence-tests) tokens-of/alphabet
  (let ((tokens "abcd"))
    (ensure (string= tokens
                     (tokens-of (make-instance 'alphabet
                                               :tokens tokens))))))

;;; Sequence constructors
(addtest (bio-sequence-tests) make-dna/rna
  (let ((seqs (list (make-dna "tagc") ; unambiguous
                    (make-dna "nnnn") ; ambiguous
                    (make-rna "uagc") ; unambiguous
                    (make-rna "nnnn"))) ; ambiguous
        (class-names (list 'dna-sequence
                           'dna-sequence
                           'rna-sequence
                           'rna-sequence))
        (alphabets (list *dna*
                         *dna*
                         *rna*
                         *rna*)))
    (mapcar #'(lambda (seq class-name alphabet)
                (ensure (subtypep (class-name (class-of seq)) class-name))
                (ensure-same alphabet (alphabet-of seq)))
            seqs class-names alphabets))
  (ensure-error
    (make-dna ""))
  (ensure-error
    (make-dna "u"))
  (ensure-error
    (make-rna "t"))
  (ensure-error
    (make-dna '(#\t #\a #\g #\c))))

(addtest (bio-sequence-tests) phred-quality
  (let ((pvals '(0.1 0.01 0.001 0.0001 0.00001))
        (quals '(10 20 30 40 50)))
    (loop
       for pval in pvals
       for qual in quals
       do (ensure-same qual (phred-quality pval)))))

(addtest (bio-sequence-tests) make-dna-quality
  (let ((seqs (list (make-dna-quality ; unambiguous
                     "agaatattctgaccccagttactttcaaga"
                     "<<<<<<<<<<<<<<<<<<<<<735513;3<"
                     :metric :phred)
                    (make-dna-quality ; ambiguous
                     "ntgccaaaaaatagaaaagtcancgatatt"
                     "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
                     :metric :phred)))
        (class-names (list 'dna-quality-sequence
                           'dna-quality-sequence)))
    (mapcar #'(lambda (seq class-name)
                (ensure (subtypep (class-name (class-of seq)) class-name)))
            seqs class-names))
  (ensure-error
   (make-dna-quality
    "agaatattctgaccccagttactttcaaga"
    "<<<<<<<"
    :metric :phred))
  (ensure-error
    (make-dna-quality
     "agaatattctgaccccagttactttcaaga"
     "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
     :metric :invalid-metric)))

;;; Utility methods
(addtest (bio-sequence-tests) simplep/string
  (ensure (simplep "acgt" (find-alphabet :dna)))
  (ensure (not (simplep "acgn" (find-alphabet :dna))))
  (ensure (simplep "acgu" (find-alphabet :rna)))
  (ensure (not (simplep "acgn" (find-alphabet :rna)))))

(addtest (bio-sequence-tests) anonymousp
  (ensure (anonymousp (make-dna "tagc")))
  (ensure (not (anonymousp (make-dna "tagc" :identity "test")))))

;;; Sequence accessors
(addtest (bio-sequence-tests) length-of/bio-sequence
  (let ((len 10))
    (ensure (= len (length-of (make-dna "aaaaaaaaaa"))))
    (ensure (= len (length-of (make-dna "aaaaaaaaaa"))))))

(addtest (bio-sequence-tests) residue-of/dna-sequence
  (let ((residues "tttt")
        (seq (make-dna "aaaa")))
    (dotimes (n (length residues))
      (setf (residue-of seq n) (aref residues n))
      (ensure (char= (residue-of seq n) (aref residues n))))))

(addtest (bio-sequence-tests) residue-of/rna-sequence
  (let ((residues "uuuu")
        (rna-seq (make-rna "aaaa")))
    (dotimes (n (length residues))
      (setf (residue-of rna-seq n) (aref residues n))
      (ensure (char= (residue-of rna-seq n) (aref residues n))))))

;;; Sequence transformations
(addtest (bio-sequence-tests) to-string/dna-sequence
  (let* ((residues "acgt")
         (seq (make-dna residues)))
    ;; no args
    (ensure (string= residues (to-string seq)))
    ;; optional arg start
     (dotimes (n 4)
       (ensure (string= (subseq residues n)
                        (to-string seq :start n))))
     (dotimes (n 4)
       (ensure (string= (subseq residues 0 n)
                        (to-string seq :start 0 :end n))))))

(addtest (bio-sequence-tests) to-string/rna-sequence
  (let* ((residues "acgu")
         (seq (make-rna residues)))
    ;; no args
    (ensure (string= residues (to-string seq)))
    ;; optional arg start
     (dotimes (n 4)
       (ensure (string= (subseq residues n)
                        (to-string seq :start n))))
     (dotimes (n 4)
       (ensure (string= (subseq residues 0 n)
                        (to-string seq :start 0 :end n))))))

(addtest (bio-sequence-tests) reverse-sequence/dna-sequence
  (let* ((residues "aaccggtt")
         (seq (make-dna residues)))
    (ensure (string= (to-string (reverse-sequence seq))
                     (reverse residues)))))

(addtest (bio-sequence-tests) reverse-sequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (seq (make-dna-quality residues quality :metric :phred))
         (rseq (reverse-sequence seq)))
    (ensure (string= (to-string rseq)
                     (reverse residues)))
    (ensure (string= (reverse quality)
                     (map-into (make-string (length residues))
                               #'encode-phred-quality
                               (quality-of rseq))))))

(addtest (bio-sequence-tests) nreverse-sequence/dna-sequence
  (let* ((residues "aaccggtt")
         (seq (make-dna residues)))
    (ensure (string= (to-string (nreverse-sequence seq))
                     (reverse residues)))))

(addtest (bio-sequence-tests) nreverse-sequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (rresidues (reverse residues))
         (rquality (reverse quality))
         (seq (make-dna-quality residues quality :metric :phred))
         (rseq (nreverse-sequence seq)))
    (ensure (string= rresidues (to-string rseq)))
    (ensure (equalp rquality
                    (map-into (make-string (length residues))
                              #'encode-phred-quality
                              (quality-of rseq))))))

(addtest (bio-sequence-tests) complement-sequence/dna-sequence
  (let ((seq (make-dna dna-residues)))
    (ensure (string= (to-string (complement-sequence seq))
                     "atcgyrmkswvhdbn-"))))

(addtest (bio-sequence-tests) complement-sequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (seq (make-dna-quality residues quality :metric :phred)))
    (ensure (string= "tcttataagactggggtcaatgaaagttct"
                     (to-string (complement-sequence seq))))
    (ensure (string= quality
                     (map-into (make-string (length residues))
                               #'encode-phred-quality
                               (quality-of seq))))))

(addtest (bio-sequence-tests) reverse-complement/dna-sequence
  (let ((seq (make-dna dna-residues)))
    (ensure (string= (to-string (reverse-complement seq))
                     "-nbdhvwskmrygcta"))))

(addtest (bio-sequence-tests) reverse-complement/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (rquality (reverse quality))
         (seq (make-dna-quality residues quality :metric :phred)))
    (ensure (string= "tcttgaaagtaactggggtcagaatattct"
                     (to-string (reverse-complement seq))))
    (let ((rcquality (quality-of (reverse-complement seq))))
      (ensure (string= rquality
                       (map-into (make-string (length residues))
                                 #'encode-phred-quality rcquality))))))

(addtest (bio-sequence-tests) nreverse-complement/dna-sequence
  (let ((seq (make-dna dna-residues)))
    (ensure (string= (to-string (nreverse-complement seq))
                     "-nbdhvwskmrygcta"))))

(addtest (bio-sequence-tests) nreverse-complement/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (rquality (reverse quality))
         (seq (make-dna-quality residues quality :metric :phred))
         (rseq (nreverse-complement seq)))
    (ensure (string= "tcttgaaagtaactggggtcagaatattct" (to-string rseq)))
    (let ((rcquality (quality-of rseq)))
      (ensure (string= rquality
                       (map-into (make-string (length residues))
                                 #'encode-phred-quality rcquality))))))

(addtest (bio-sequence-tests) subsequence/bio-sequence
  (let* ((residues "aaggccttaaggcctt")
         (seq (make-dna residues)))
    (ensure (string= (subseq residues 0 5) ; aaggc
                     (to-string (subsequence seq 0 5))))
    (ensure (string= residues
                     (to-string (subsequence seq 0))))))

(addtest (bio-sequence-tests) subsequence/dna-quality-sequence
  (let* ((residues "agaatattctgaccccagttactttcaaga")
         (quality "<<<<<<<<<<<<<<<<<<<<<735513;3<")
         (seq (make-dna-quality residues quality :metric :phred)))
    (ensure (string= (subseq residues 0 5)
                     (to-string (subsequence seq 0 5))))
    (loop
       for q across (quality-of (subsequence seq 0 5))
       do (ensure (= 27 q)))))

(addtest (bio-sequence-tests) residue-frequencies/bio-sequence
  (let ((seq (make-dna "aacccgggt")))
    (mapc #'(lambda (token freq)
              (ensure (eq (cdr (assoc token (residue-frequencies
                                             seq
                                             (find-alphabet :dna))
                                      :test #'char= )) freq)))
          '(#\a #\c #\g #\t)
          '(2 3 3 1))))

(addtest (bio-sequence-tests) search-sequence
  (let ((seq (make-dna "tacgagtcgttttagcgcgattatataa"))
        (sub (make-dna "gtttt")))
    (ensure (= 8 (search-sequence sub seq)))
    (ensure (= 9 (search-sequence sub seq :start1 1)))
    (ensure (= 8 (search-sequence sub seq :start2 1)))
    (ensure (= 3 (search-sequence sub seq :end1 1)))
    (ensure-null (search-sequence sub seq :end2 1))))
