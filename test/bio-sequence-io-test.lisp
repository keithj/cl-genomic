;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
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
  (with-seq-input (seqi filespec format :alphabet :dna :metric :sanger)
    (loop
       for seq = (next seqi)
       count seq into total
       while (has-more-p seqi)
       finally (return total))))

(defmacro with-test-file ((stream filespec) &body body)
  (with-gensyms (fs)
    `(let ((,fs (merge-pathnames ,filespec)))
      (with-seqi (,stream ,fs)
        ,@body))))

(defmacro with-test-mapped-seq ((mseq seq filespec) &body body)
  (with-gensyms (tmp-filespec seqi)
    `(let* ((,seq (with-seq-input (,seqi ,filespec :fasta)
                    (next ,seqi)))
            (,tmp-filespec (tmp-pathname :tmpdir (merge-pathnames "data"))))
      (write-pure-sequence ,seq ,tmp-filespec)
      (with-mapped-dna (,mseq :filespec ,tmp-filespec
                              :length (length-of ,seq))
        ,@body)
      (delete-file ,tmp-filespec))))


(deftestsuite bio-sequence-io-tests (cl-genomic-tests)
  ())

(addtest (bio-sequence-io-tests) fasta/1
  (with-seq-input (seqi (merge-pathnames "data/simple-dna1.fasta")
                        :fasta :alphabet :dna)
    (let ((seq (next seqi)))
      (ensure (subtypep (type-of seq) 'dna-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (not (virtualp seq)))
      (ensure (= 210 (length-of seq)))
      (ensure (string= "Test1" (identity-of seq))))))

(addtest (bio-sequence-io-tests) fasta/2
  (with-seq-input (seqi (merge-pathnames "data/simple-dna1.fasta")
                        :fasta :alphabet :dna :virtual t)
    (let ((seq (next seqi)))
      (ensure (subtypep (type-of seq) 'dna-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (virtualp seq))
      (ensure (= 210 (length-of seq)))
      (ensure (string= "Test1" (identity-of seq))))))

(addtest (bio-sequence-io-tests) fasta/3
  (with-seq-input (seqi (merge-pathnames "data/iupac-dna1.fasta")
                        :fasta :alphabet :dna)
    (let ((seq (next seqi)))
      (ensure (subtypep (type-of seq) 'dna-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (= 210 (length-of seq)))
      (ensure (string= "Test1" (identity-of seq))))))

(addtest (bio-sequence-io-tests) fasta/4
  (with-seq-input (seqi (merge-pathnames "data/simple-dna2.fasta")
                        :fasta :alphabet :dna)
    (dotimes (n 2)
      (let ((seq (next seqi)))
        (ensure (subtypep (type-of seq) 'dna-sequence))
        (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
        (ensure (= 280 (length-of seq)))
        (ensure (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
    (ensure-null (next seqi))))

(addtest (bio-sequence-io-tests) fasta/5
  (with-seq-input (seqi (merge-pathnames "data/simple-dna2.fasta")
                        :fasta :alphabet :dna :virtual t)
    (dotimes (n 2)
      (let ((seq (next seqi)))
        (ensure (subtypep (type-of seq) 'dna-sequence))
        (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
        (ensure (virtualp seq))
        (ensure (= 280 (length-of seq)))
        (ensure (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
    (ensure-null (next seqi))))

(addtest (bio-sequence-io-tests) fasta/6
  (with-seq-input (seqi (merge-pathnames "data/iupac-dna2.fasta")
                        :fasta :alphabet :dna)
    (dotimes (n 2)
      (let ((seq (next seqi)))
        (ensure (subtypep (type-of seq) 'dna-sequence))
        (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
        (ensure (= 280 (length-of seq)))
        (ensure (string= (format nil "Test~a" (1+ n)) (identity-of seq)))))
    (ensure-null (next seqi))))

(addtest (bio-sequence-io-tests) fasta/7
  (with-seq-input (seqi (merge-pathnames "data/iupac-dna2.fasta")
                        :fasta :alphabet :dna :parser (make-instance
                                                       'raw-sequence-parser))
    (dotimes (n 2)
      (let ((seq (next seqi)))
        (ensure (eql :dna (assocdr :alphabet seq)))
        (ensure (= 280 (length (assocdr :residues seq))))
        (ensure (string= (format nil "Test~a" (1+ n))
                         (assocdr :identity seq)))))
    (ensure-null (next seqi))))

(addtest (bio-sequence-io-tests) fasta/8
  (with-seq-input (seqi (merge-pathnames "data/phred.fastq")
                        :fasta :alphabet :dna) ; fastq!
    (ensure-condition malformed-record-error
      (next seqi))))

(addtest (bio-sequence-io-tests) write-fasta-sequence/1
  (let ((seq (make-dna "acgtn" :identity "foo"))
        (tmp-filespec (tmp-pathname :tmpdir (merge-pathnames "data")
                                    :type "fa")))
    (dolist (args '((nil "acgtn")
                    (:lower "acgtn")
                    (:upper "ACGTN")))
      (with-open-file (stream tmp-filespec :direction :io
                              :element-type 'base-char
                              :external-format :ascii)
        (write-fasta-sequence seq stream :token-case (first args))
        (finish-output stream)
        (file-position stream 0)
        (ensure (string= ">foo" (read-line stream)))
        (ensure (string= (second args) (read-line stream))))
      (delete-file tmp-filespec))))

(addtest (bio-sequence-io-tests) fastq/1
  (with-seq-input (seqi (merge-pathnames "data/phred.fastq")
                        :fastq :alphabet :dna)
    (do ((seq (next seqi) (next seqi)))
        ((null seq) t)
      (ensure (subtypep (type-of seq) 'dna-quality-sequence))
      (ensure (eql (find-alphabet :dna) (alphabet-of seq)))
      (ensure (= 35 (length-of seq)))
      (ensure (string= "IL13" (identity-of seq) :start2 0 :end2 4)))))

(addtest (bio-sequence-io-tests) fastq/2
  (with-seq-input (seqi (merge-pathnames "data/phred.fastq")
                        :fastq :alphabet :dna :parser (make-instance
                                                       'raw-sequence-parser))
    (do ((seq (next seqi) (next seqi)))
        ((null seq) t)
      (ensure (eql :dna (assocdr :alphabet seq)))
      (ensure (= 35 (length (assocdr :residues seq))))
      (ensure (= 35 (length (assocdr :quality seq))))
      (ensure (string= "IL13" (assocdr :identity seq) :start2 0 :end2 4)))))

(addtest (bio-sequence-io-tests) fastq/3
  (with-seq-input (seqi (merge-pathnames "data/simple-dna1.fasta") ; fasta!
                        :fastq :alphabet :dna :metric :sanger)
    (ensure-condition malformed-record-error
      (next seqi))))

(addtest (bio-sequence-io-tests) fastq/4
  ;; Identifier is an empty string. This is technically legal because
  ;; the Fastq paper explicitly states that there is no length limit
  ;; on the title field. The authors probably meant no upper length
  ;; limit only.
  (with-seq-input (seqi (merge-pathnames "data/no_identifier.fastq")
                        :fastq :alphabet :dna :metric :sanger)
    (ensure (string= "" (identity-of (next seqi))))
    (let ((seq2 (next seqi)))
      (ensure (string= "" (identity-of seq2)))
      ;; Currently we ignore the description for Fastq
      ;; (ensure (string= "a description" (description-of seq2)))
      (ensure (null (description-of seq2))))))

(addtest (bio-sequence-io-tests) write-fastq-sequence/1
  (let ((seq (make-dna-quality "acgtn" "<<<<<"
                               :identity "foo"
                               :metric :sanger))
        (tmp-filespec (tmp-pathname :tmpdir (merge-pathnames "data")
                                    :type "fq")))
    (dolist (args '((nil "acgtn")
                    (:lower"acgtn")
                    (:upper "ACGTN")))
      (with-open-file (stream tmp-filespec :direction :io
                              :element-type 'base-char
                              :external-format :ascii)
        (write-fastq-sequence seq stream :token-case (first args))
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
                         (pathname-extender filespec
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
                         (pathname-extender filespec
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
        (out-filespec (namestring
                       (tmp-pathname :tmpdir (merge-pathnames "data")
                                     :type "fasta"))))
    (convert-sequence-file in-filespec :fastq out-filespec :fasta)
    (with-seq-input (fqi in-filespec :fastq :alphabet :dna :metric :sanger)
      (with-seq-input (fai out-filespec :fasta :alphabet :dna)
        (ensure (loop
                   for fq = (next fqi)
                   while fq
                   always (let ((fa (next fai)))
                            (and (string= (identity-of fq)
                                          (identity-of fa))
                                 (string= (coerce-sequence fq 'string)
                                          (coerce-sequence fa 'string))))))
        (ensure-null (next fai))))
    (delete-file out-filespec)))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/1
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (dna-sequence-p mseq))
    (ensure (double-stranded-p mseq))
    (ensure (simplep mseq))
    (ensure (not (ambiguousp mseq)))
    (ensure (= (length-of seq) (length-of mseq)))
     (dotimes (n (length-of seq))
       (setf (residue-of mseq n) #\n)
       (ensure (char= #\n (residue-of mseq n))))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/2
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (let ((str (coerce-sequence seq 'string)))
    ;; no args
    (ensure (string= str (coerce-sequence mseq 'string)))
    ;; optional arg start
    (dotimes (n (length str))
      (ensure (string= (subseq str n)
                       (coerce-sequence mseq 'string :start n))))
    ;; optional arg start end
    (dotimes (n (length str))
      (ensure (string= (subseq str 0 n)
                       (coerce-sequence mseq 'string :start 0 :end n)))))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/3
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (string= (coerce-sequence (reverse-sequence seq) 'string)
                     (coerce-sequence (reverse-sequence mseq) 'string)))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/4
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (string= (coerce-sequence (nreverse-sequence seq) 'string)
                     (coerce-sequence (nreverse-sequence mseq) 'string)))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/5
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (string= (coerce-sequence (complement-sequence seq) 'string)
                     (coerce-sequence (complement-sequence mseq) 'string)))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/6
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (string= (coerce-sequence (ncomplement-sequence seq) 'string)
                     (coerce-sequence (ncomplement-sequence mseq) 'string)))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/7
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (string= (coerce-sequence (reverse-complement seq) 'string)
                     (coerce-sequence (reverse-complement mseq) 'string)))))

(addtest (bio-sequence-io-tests) mapped-dna-sequence/8
  (with-test-mapped-seq (mseq seq (merge-pathnames "data/simple-dna1.fasta"))
    (ensure (string= (coerce-sequence (nreverse-complement seq) 'string)
                     (coerce-sequence (nreverse-complement mseq) 'string)))))
