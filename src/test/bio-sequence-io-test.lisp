;;;
;;; Copyright (C) 2007-2009 Keith James. All rights reserved.
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

(defun count-seq-records (filespec format)
  (with-ascii-li-stream (stream filespec)
    (let ((gen (make-seq-input stream format :alphabet :dna
                                             :metric :phred)))
      (loop
         as seq = (next gen)
         count 1 into total
         while (has-more-p gen)
         finally (return total)))))

(defmacro with-test-file ((stream filespec) &body body)
  (with-gensyms (fs)
    `(let ((,fs (merge-pathnames ,filespec)))
      (with-ascii-li-stream (,stream ,fs)
        ,@body))))

(deftestsuite bio-sequence-io-tests (cl-genomic-tests)
  ())

(addtest (bio-sequence-io-tests) fasta/1
  (with-test-file (stream "data/simple-dna1.fasta")
    (let* ((gen (make-seq-input stream :fasta :alphabet :dna))
           (cur (current gen))
           (seq (next gen)))
      (ensure (eql cur seq))
      (ensure (subtypep (type-of seq) 'dna-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (not (virtualp seq)))
      (ensure (= 210 (length-of seq)))
      (ensure (string= "Test1" (identity-of seq))))))

(addtest (bio-sequence-io-tests) fasta/2
  (with-test-file (stream "data/simple-dna1.fasta")
    (let* ((gen (make-seq-input stream :fasta :alphabet :dna :virtual t))
           (cur (current gen))
           (seq (next gen)))
      (ensure (eql cur seq))
      (ensure (subtypep (type-of seq) 'dna-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (virtualp seq))
      (ensure (= 210 (length-of seq)))
      (ensure (string= "Test1" (identity-of seq))))))

(addtest (bio-sequence-io-tests) fasta/3
  (with-test-file (stream "data/iupac-dna1.fasta")
    (let* ((gen (make-seq-input stream :fasta :alphabet :dna))
           (seq (next gen)))
      (ensure (subtypep (type-of seq) 'dna-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (= 210 (length-of seq)))
      (ensure (string= "Test1" (identity-of seq))))))

(addtest (bio-sequence-io-tests) fasta/4
  (with-test-file (stream "data/simple-dna2.fasta")
    (let ((gen (make-seq-input stream :fasta :alphabet :dna)))
      (dotimes (n 2)
        (let ((cur (current gen))
              (seq (next gen)))
          (ensure (eql cur seq))
          (ensure (subtypep (type-of seq) 'dna-sequence))
          (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
          (ensure (= 280 (length-of seq)))
          (ensure (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
      (ensure-null (next gen)))))

(addtest (bio-sequence-io-tests) fasta/5
  (with-test-file (stream "data/simple-dna2.fasta")
    (let ((gen (make-seq-input stream :fasta :alphabet :dna :virtual t)))
      (dotimes (n 2)
        (let ((cur (current gen))
              (seq (next gen)))
          (ensure (eql cur seq))
          (ensure (subtypep (type-of seq) 'dna-sequence))
          (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
          (ensure (virtualp seq))
          (ensure (= 280 (length-of seq)))
          (ensure (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
      (ensure-null (next gen)))))

(addtest (bio-sequence-io-tests) fasta/6
  (with-test-file (stream "data/iupac-dna2.fasta")
    (let ((gen (make-seq-input stream :fasta :alphabet :dna)))
      (dotimes (n 2)
        (let ((cur (current gen))
              (seq (next gen)))
          (ensure (eql cur seq))
          (ensure (subtypep (type-of seq) 'dna-sequence))
          (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
          (ensure (= 280 (length-of seq)))
          (ensure (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
      (ensure-null (next gen)))))

(addtest (bio-sequence-io-tests) fasta/7
  (with-test-file (stream "data/iupac-dna2.fasta")
    (let ((gen (make-seq-input stream :fasta :alphabet :dna
                                      :parser (make-instance
                                               'raw-sequence-parser))))
      (dotimes (n 2)
        (let ((seq (next gen)))
          (ensure (eql :dna (gpu:assocdr :alphabet seq)))
          (ensure (= 280 (length (gpu:assocdr :residues seq))))
          (ensure (string= (format nil "Test~a" (1+ n))
                           (gpu:assocdr :identity seq)))))
      (ensure-null (next gen)))))

(addtest (bio-sequence-io-tests) fasta/8
  (with-test-file (stream "data/phred.fastq") ; fastq!
    (ensure-condition malformed-record-error
        (make-seq-input (make-line-input-stream stream) :fasta
                        :alphabet :dna))))

(addtest (bio-sequence-io-tests) write-fasta-sequence/1
  (let ((seq (make-dna "acgtn" :identity "foo"))
        (tmp-filespec (iou:make-tmp-pathname
                       :tmpdir (merge-pathnames "data")
                       :type "fa")))
    (dolist (args '((nil "acgtn")
                    (:lower "acgtn")
                    (:upper "ACGTN")))
      (with-open-file (stream tmp-filespec :direction :io
                              :element-type 'base-char
                              :external-format :ascii)
        (bs::write-fasta-sequence seq stream :token-case (first args))
        (finish-output stream)
        (file-position stream 0)
        (ensure (string= ">foo" (read-line stream)))
        (ensure (string= (second args) (read-line stream))))
      (delete-file tmp-filespec))))

(addtest (bio-sequence-io-tests) fastq/1
  (with-test-file (stream "data/phred.fastq")
    (let ((gen (make-seq-input stream :fastq :alphabet :dna)))
      (do ((cur (current gen) (current gen))
           (seq (next gen) (next gen)))
          ((null seq) t)
        (ensure (eql cur seq))
        (ensure (subtypep (type-of seq) 'dna-quality-sequence))
        (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
        (ensure (= 35 (length-of seq)))
        (ensure (string= "IL13" (identity-of seq) :start2 0 :end2 4))))))

(addtest (bio-sequence-io-tests) fastq/2
  (with-test-file (stream "data/phred.fastq")
    (let ((gen (make-seq-input stream :fastq :alphabet :dna
                                      :parser (make-instance
                                               'raw-sequence-parser))))
       (do ((cur (current gen) (current gen))
            (seq (next gen) (next gen)))
           ((null seq) t)
         (ensure (eql cur seq))
         (ensure (eql :dna (gpu:assocdr :alphabet seq)))
         (ensure (= 35 (length (gpu:assocdr :residues seq))))
         (ensure (= 35 (length (gpu:assocdr :quality seq))))
         (ensure (string= "IL13" (gpu:assocdr :identity seq)
                          :start2 0 :end2 4))))))

(addtest (bio-sequence-io-tests) fastq/3
  (with-test-file (stream "data/simple-dna1.fasta") ; fasta!
    (ensure-condition malformed-record-error
      (make-seq-input (make-line-input-stream stream) :fastq
                      :alphabet :dna
                      :metric :phred))))

(addtest (bio-sequence-io-tests) write-fastq-sequence/1
  (let ((seq (make-dna-quality "acgtn" "<<<<<"
                               :identity "foo"
                               :metric :phred))
        (tmp-filespec (iou:make-tmp-pathname
                       :tmpdir (merge-pathnames "data")
                       :type "fq")))
    (dolist (args '((nil "acgtn")
                    (:lower"acgtn")
                    (:upper "ACGTN")))
      (with-open-file (stream tmp-filespec :direction :io
                              :element-type 'base-char
                              :external-format :ascii)
        (bs::write-fastq-sequence seq stream :token-case (first args))
        (finish-output stream)
        (file-position stream 0)
        (ensure (string= "@foo" (read-line stream)))
        (ensure (string= (second args) (read-line stream)))
        (ensure (string= "+" (read-line stream)))
        (ensure (string= "<<<<<" (read-line stream))))
      (delete-file tmp-filespec))))

(addtest (bio-sequence-io-tests) split-sequence-file/1
  (let ((filespec (namestring (merge-pathnames "data/split-test-dna1.fasta"))))
    (split-sequence-file filespec :fasta
                         (make-pathname-ext filespec
                                            :type "fasta" :separator #\.)
                         :chunk-size 2))
  (dolist (args '(("data/split-test-dna1.0.fasta" 2)
                  ("data/split-test-dna1.1.fasta" 2)
                  ("data/split-test-dna1.2.fasta" 2)
                  ("data/split-test-dna1.3.fasta" 1)))
    (let ((chunk (merge-pathnames (first args))))
      (ensure (= (second args) (count-seq-records chunk :fasta)))
      (delete-file chunk))))

(addtest (bio-sequence-io-tests) split-sequence-file/2
  (let ((filespec (namestring (merge-pathnames "data/phred.fastq"))))
    (split-sequence-file filespec :fastq
                         (make-pathname-ext filespec
                                            :type "fastq" :separator #\.)
                         :chunk-size 2))
  (dolist (args '(("data/phred.0.fastq" 2)
                  ("data/phred.1.fastq" 2)
                  ("data/phred.2.fastq" 2)
                  ("data/phred.3.fastq" 1)))
    (let ((chunk (merge-pathnames (first args))))
      (ensure (= (second args) (count-seq-records chunk :fastq)))
      (delete-file chunk))))

(addtest (bio-sequence-io-tests) convert-sequence-file/1
  (let ((in-filespec (merge-pathnames "data/phred.fastq"))
        (out-filespec (namestring (iou:make-tmp-pathname
                                   :tmpdir (merge-pathnames "data")
                                   :type "fasta"))))
    (convert-sequence-file in-filespec :fastq out-filespec :fasta)
    (with-ascii-li-stream (fqs in-filespec)
      (with-ascii-li-stream (fas out-filespec)
        (let ((fq-gen (make-seq-input fqs :fastq :alphabet :dna
                                                 :metric :phred))
              (fa-gen (make-seq-input fas :fasta :alphabet :dna)))
          (ensure (loop
                     as fq = (next fq-gen)
                     while fq
                     always (let ((fa (next fa-gen)))
                              (and (string= (identity-of fq)
                                            (identity-of fa))
                                   (string= (coerce-sequence fq 'string)
                                            (coerce-sequence fa 'string))))))
          (ensure-null (next fa-gen)))))
    (delete-file out-filespec)))
