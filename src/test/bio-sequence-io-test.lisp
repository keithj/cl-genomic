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

(in-package :cl-bio-test)

(fiveam:in-suite cl-bio-system:testsuite)

;;; Test reading unambiguous/IUPAC Fasta DNA/RNA
(test read-bio-sequence/fasta/dna-simple
  (with-open-file (fs (merge-pathnames "data/simple-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((stream (make-line-input-stream fs))
           (s (read-bio-sequence stream :fasta
                                 :alphabet :dna
                                 :ambiguity nil)))
      (is (eql 'simple-dna-sequence (type-of s)))
      (is (eql (find-alphabet :dna) (alphabet-of s)))
      (is-false (virtualp s))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (identity-of s))))))

(test read-bio-sequence/fasta/dna-simple/virtual
  (with-open-file (fs (merge-pathnames "data/simple-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((stream (make-line-input-stream fs))
           (s (read-bio-sequence stream :fasta
                                 :alphabet :dna
                                 :ambiguity nil
                                 :virtualp t)))
      (is (eql 'simple-dna-sequence (type-of s)))
      (is (eql (find-alphabet :dna) (alphabet-of s)))
      (is-true (virtualp s))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (identity-of s))))))

(test read-bio-sequence/fasta/dna-iupac
  (with-open-file (fs (merge-pathnames "data/iupac-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((stream (make-line-input-stream fs))
           (s (read-bio-sequence stream :fasta
                                 :alphabet :dna
                                 :ambiguity :iupac)))
      (is (eql 'iupac-dna-sequence (type-of s)))
      (is (eql (find-alphabet :iupac-dna) (alphabet-of s)))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (identity-of s))))))

(test read-bio-sequence/fasta/dna-auto
  (with-open-file (fs (merge-pathnames "data/iupac-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((stream (make-line-input-stream fs))
           (s (read-bio-sequence stream :fasta
                                 :alphabet :dna
                                 :ambiguity :auto)))
      (is (eql 'iupac-dna-sequence (type-of s)))
      (is (eql (find-alphabet :iupac-dna) (alphabet-of s)))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (identity-of s))))))

(test read-bio-sequence/multifasta/dna-simple
  (with-open-file (fs (merge-pathnames "data/simple-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((stream (make-line-input-stream fs)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence stream :fasta
                                    :alphabet :dna
                                    :ambiguity nil)))
          (is (eql 'simple-dna-sequence (type-of s)))
          (is (eql (find-alphabet :dna) (alphabet-of s)))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of s)))))
      (is (null (read-bio-sequence stream :fasta
                                   :alphabet :dna))))))

(test read-bio-sequence/multifasta/dna-simple/virtual
  (with-open-file (fs (merge-pathnames "data/simple-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((stream (make-line-input-stream fs)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence stream :fasta
                                    :alphabet :dna
                                    :ambiguity nil
                                    :virtualp t)))
          (is (eql 'simple-dna-sequence (type-of s)))
          (is (eql (find-alphabet :dna) (alphabet-of s)))
          (is-true (virtualp s))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of s)))))
      (is (null (read-bio-sequence stream :fasta
                                   :alphabet :dna))))))

(test read-bio-sequence/multifasta/dna-iupac
  (with-open-file (fs (merge-pathnames "data/iupac-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((stream (make-line-input-stream fs)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence stream :fasta
                                    :alphabet :dna
                                    :ambiguity :iupac)))
          (is (eql 'iupac-dna-sequence (type-of s)))
          (is (eql (find-alphabet :iupac-dna) (alphabet-of s)))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of s)))))
      (is (null (read-bio-sequence stream :fasta
                                   :alphabet :dna))))))

(test read-bio-sequence/fastq/simple
  (with-open-file (fs (merge-pathnames "data/phred.fq")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((stream (make-line-input-stream fs)))
      (do ((s (read-bio-sequence stream :fastq
                                 :alphabet :dna
                                 :ambiguity :auto
                                 :metric :phred)
              (read-bio-sequence stream :fastq
                                 :alphabet :dna
                                 :ambiguity :auto
                                 :metric :phred)))
          ((null s) t)
        (is (eql 'simple-dna-quality-sequence (type-of s)))
        (is (eql (find-alphabet :dna) (alphabet-of s)))
        (is (= 35 (length-of s)))
        (is (string= "IL13" (identity-of s) :start2 0 :end2 4))))))
