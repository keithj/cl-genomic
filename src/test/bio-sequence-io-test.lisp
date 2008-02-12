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
(test read-bio-sequence/fasta/dna
  (with-open-file (fs (merge-pathnames "data/simple-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((stream (make-line-input-stream fs))
           (s (read-bio-sequence stream
                                 :alphabet :dna
                                 :format :fasta)))
      (is (eql 'simple-dna-sequence (type-of s)))
      (is (eql *dna* (alphabet-of s)))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (identity-of s)))))
  (with-open-file (fs (merge-pathnames "data/iupac-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((stream (make-line-input-stream fs))
           (s (read-bio-sequence stream
                                 :alphabet :dna
                                 :ambiguity :iupac
                                 :format :fasta)))
      (is (eql 'iupac-dna-sequence (type-of s)))
      (is (eql *iupac-dna* (alphabet-of s)))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (identity-of s))))))

(test read-bio-sequence/multifasta/dna
  (with-open-file (fs (merge-pathnames "data/simple-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((stream (make-line-input-stream fs)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence stream
                                    :alphabet :dna
                                    :format :fasta)))
          (is (eql 'simple-dna-sequence (type-of s)))
          (is (eql *dna* (alphabet-of s)))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of s)))))
      (is (null (read-bio-sequence stream
                                   :alphabet :dna
                                   :format :fasta)))))
  (with-open-file (fs (merge-pathnames "data/iupac-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((stream (make-line-input-stream fs)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence stream
                                    :alphabet :dna
                                    :ambiguity :iupac
                                    :format :fasta)))
          (is (eql 'iupac-dna-sequence (type-of s)))
          (is (eql *iupac-dna* (alphabet-of s)))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of s)))))
      (is (null (read-bio-sequence stream
                                   :alphabet :dna
                                   :format :fasta))))))
