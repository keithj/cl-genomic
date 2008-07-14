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

(defun count-seq-records (filespec format)
  (with-open-file (stream filespec :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let ((gen (make-seq-input (make-line-input-stream stream) format
                               :alphabet :dna
                               :metric :phred)))
      (loop
         as seq = (next gen)
         count 1 into total
         while (has-more-p gen)
         finally (return total)))))


(fiveam:in-suite cl-bio-system:testsuite)

(test bio-sequence-io/fasta/dna-simple
  (with-open-file (fs (merge-pathnames "data/simple-dna1.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fasta :alphabet :dna))
           (cur (current gen))
           (seq (next gen)))
      (is (eql cur seq))
      (is (subtypep (type-of seq)  'dna-sequence))
      (is (eql (find-alphabet :dna) (alphabet-of seq)))
      (is-false (virtualp seq))
      (is (= 210 (length-of seq)))
      (is (string= "Test1" (identity-of seq))))))

(test bio-sequence-io/fasta/dna-simple/virtual
  (with-open-file (fs (merge-pathnames "data/simple-dna1.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fasta :alphabet :dna :virtual t))
           (cur (current gen))
           (seq (next gen)))
      (is (eql cur seq))
      (is (subtypep (type-of seq) 'dna-sequence))
      (is (eql (find-alphabet :dna) (alphabet-of seq)))
      (is-true (virtualp seq))
      (is (= 210 (length-of seq)))
      (is (string= "Test1" (identity-of seq))))))

(test bio-sequence-io/fasta/dna-iupac
  (with-open-file (fs (merge-pathnames "data/iupac-dna1.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fasta :alphabet :dna))
           (seq (next gen)))
      (is (subtypep (type-of seq) 'dna-sequence))
      (is (eql (find-alphabet :dna) (alphabet-of seq)))
      (is (= 210 (length-of seq)))
      (is (string= "Test1" (identity-of seq))))))

(test bio-sequence-io/multifasta/dna-simple
  (with-open-file (fs (merge-pathnames "data/simple-dna2.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fasta :alphabet :dna)))
      (dotimes (n 2)
        (let ((cur (current gen))
              (seq (next gen)))
          (is (eql cur seq))
          (is (subtypep (type-of seq) 'dna-sequence))
          (is (eql (find-alphabet :dna) (alphabet-of seq)))
          (is (= 280 (length-of seq)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
      (is (null (next gen))))))

(test bio-sequence-io/multifasta/dna-simple/virtual
  (with-open-file (fs (merge-pathnames "data/simple-dna2.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fasta :alphabet :dna :virtual t)))
      (dotimes (n 2)
        (let ((cur (current gen))
              (seq (next gen)))
          (is (eql cur seq))
          (is (subtypep (type-of seq) 'dna-sequence))
          (is (eql (find-alphabet :dna) (alphabet-of seq)))
          (is-true (virtualp seq))
          (is (= 280 (length-of seq)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
      (is (null (next gen))))))

(test bio-sequence-io/multifasta/dna-iupac
  (with-open-file (fs (merge-pathnames "data/iupac-dna2.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fasta :alphabet :dna)))
      (dotimes (n 2)
        (let ((cur (current gen))
              (seq (next gen)))
          (is (eql cur seq))
          (is (subtypep (type-of seq) 'dna-sequence))
          (is (eql (find-alphabet :dna) (alphabet-of seq)))
          (is (= 280 (length-of seq)))
          (is (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
      (is (null (next gen))))))

(test bio-sequence-io/multifasta/raw
  (with-open-file (fs (merge-pathnames "data/iupac-dna2.fa")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (parser (make-instance 'raw-sequence-parser))
           (gen (make-seq-input stream :fasta :alphabet :dna
                                :parser parser)))
      (dotimes (n 2)
        (let ((seq (next gen)))
          (is (eql :dna (gpu:assocdr :alphabet seq)))
          (is (= 280 (length (gpu:assocdr :residues seq))))
          (is (string= (format nil "Test~a" (1+ n))
                       (gpu:assocdr :identity seq)))))
      (is (null (next gen))))))

(test bio-sequence-io/fastq/simple
  (with-open-file (fs (merge-pathnames "data/phred.fq")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (gen (make-seq-input stream :fastq :alphabet :dna)))
      (do ((cur (current gen) (current gen))
           (seq (next gen) (next gen)))
          ((null seq) t)
        (is (eql cur seq))
        (is (subtypep (type-of seq) 'dna-quality-sequence))
        (is (eql (find-alphabet :dna) (alphabet-of seq)))
        (is (= 35 (length-of seq)))
        (is (string= "IL13" (identity-of seq) :start2 0 :end2 4))))))

(test bio-sequence-io/fastq/raw
  (with-open-file (fs (merge-pathnames "data/phred.fq")
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let* ((stream (make-line-input-stream fs))
           (parser (make-instance 'raw-sequence-parser))
           (gen (make-seq-input stream :fastq :alphabet :dna
                                :parser parser)))
       (do ((cur (current gen) (current gen))
            (seq (next gen) (next gen)))
           ((null seq) t)
         (is (eql cur seq))
         (is (eql :dna (gpu:assocdr :alphabet seq)))
         (is (= 35 (length (gpu:assocdr :residues seq))))
         (is (= 35 (length (gpu:assocdr :quality seq))))
         (is (string= "IL13" (gpu:assocdr :identity seq)
                       :start2 0 :end2 4))))))

(test bio-sequence-io/fasta/error
  (with-open-file (fs (merge-pathnames "data/phred.fq") ; fastq!
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (signals bio-sequence-io-error
        (make-seq-input (make-line-input-stream fs) :fasta
                        :alphabet :dna))))

(test bio-sequence-io/fastq/error
  (with-open-file (fs (merge-pathnames "data/simple-dna1.fa") ; fasta!
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (signals bio-sequence-io-error
      (make-seq-input (make-line-input-stream fs) :fastq
                      :alphabet :dna
                      :metric :phred))))

(test write-fasta-sequence
  (let ((seq (make-dna "acgtn" :identity "foo"))
        (tmp-filespec (iou:make-tmp-pathname
                       :tmpdir (merge-pathnames "data")
                       :type "fa")))
    (dolist (args '((nil "acgtn") ( :lowercase "acgtn") (:uppercase "ACGTN")))
      (with-open-file (stream tmp-filespec :direction :io
                       :element-type 'base-char
                       :external-format :ascii)
        (bio-sequence::write-fasta-sequence seq stream
                                            :token-case (first args))
        (finish-output stream)
        (file-position stream 0)
        (is (string= ">foo" (read-line stream)))
        (is (string= (second args) (read-line stream))))
      (delete-file tmp-filespec))))

(test write-fastq-sequence
  (let ((seq (make-dna-quality "acgtn" "<<<<<"
                               :identity "foo"
                               :metric :phred))
        (tmp-filespec (iou:make-tmp-pathname
                       :tmpdir (merge-pathnames "data")
                       :type "fq")))
    (dolist (args '((nil "acgtn") (:lowercase "acgtn") (:uppercase "ACGTN")))
      (with-open-file (stream tmp-filespec :direction :io
                       :element-type 'base-char
                       :external-format :ascii)
        (bio-sequence::write-fastq-sequence seq stream
                                            :token-case (first args))
        (finish-output stream)
        (file-position stream 0)
        (is (string= "@foo" (read-line stream)))
        (is (string= (second args) (read-line stream)))
        (is (string= "+" (read-line stream)))
        (is (string= "<<<<<" (read-line stream))))
      (delete-file tmp-filespec))))

(test split-sequence-file/fastq
  (let ((filespec (namestring (merge-pathnames "data/phred.fq"))))
    (split-sequence-file filespec :fastq
                         (make-pathname-ext filespec
                                            :type "fq" :separator #\.)
                         :chunk-size 2))
  (dolist (args '(("data/phred.0.fq" 2)
                  ("data/phred.1.fq" 2)
                  ("data/phred.2.fq" 2)
                  ("data/phred.3.fq" 1)))
    (let ((chunk (merge-pathnames (first args))))
      (is (= (second args) (count-seq-records chunk :fastq)))
      (delete-file chunk))))

(test split-sequence-file/fasta
  (let ((filespec (namestring (merge-pathnames "data/split-test-dna1.fa"))))
    (split-sequence-file filespec :fasta
                         (make-pathname-ext filespec
                                            :type "fa" :separator #\.)
                         :chunk-size 2))
  (dolist (args '(("data/split-test-dna1.0.fa" 2)
                  ("data/split-test-dna1.1.fa" 2)
                  ("data/split-test-dna1.2.fa" 2)
                  ("data/split-test-dna1.3.fa" 1)))
    (let ((chunk (merge-pathnames (first args))))
      (is (= (second args) (count-seq-records chunk :fasta)))
      (delete-file chunk))))

(test convert-sequence-file/fastq->fasta
  (let ((in-filespec (merge-pathnames "data/phred.fq"))
        (out-filespec (namestring (iou:make-tmp-pathname
                                   :tmpdir (merge-pathnames "data")
                                   :type "fa"))))
    (convert-sequence-file in-filespec :fastq out-filespec :fasta)
    (with-open-file (fqs in-filespec :direction :input)
      (with-open-file (fas out-filespec :direction :input)
        (let ((fq-gen (make-seq-input (make-line-input-stream fqs) :fastq
                                      :alphabet :dna
                                      :metric :phred))
              (fa-gen (make-seq-input (make-line-input-stream fas) :fasta
                                      :alphabet :dna)))
          (is-true (loop
                      as fq = (next fq-gen)
                      while fq
                      always (let ((fa (next fa-gen)))
                               (and (string= (identity-of fq)
                                             (identity-of fa))
                                    (string= (to-string fq)
                                             (to-string fa))))))
          (is (null (next fa-gen))))))
    (delete-file out-filespec)))
