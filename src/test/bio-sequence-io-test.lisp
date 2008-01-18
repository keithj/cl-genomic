
(in-package :cl-bio-test)

(fiveam:in-suite cl-bio-system:testsuite)

;; Test reading unambiguous/IUPAC Fasta DNA/RNA

(test read-bio-sequence/fasta/dna
  (with-open-file (stream (merge-pathnames "data/simple-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((line-buffer (make-line-buffer stream))
           (s (read-bio-sequence line-buffer
                                 :alphabet :dna
                                 :format :fasta
                                 :callback #'make-seq-from-alist)))
      (is (eql 'simple-dna-sequence (type-of s)))
      (is (eql *dna* (alphabet-of s)))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (name-of s)))))
  (with-open-file (stream (merge-pathnames "data/iupac-dna1.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((line-buffer (make-line-buffer stream))
           (s (read-bio-sequence line-buffer
                                 :alphabet :dna
                                 :ambiguity :iupac
                                 :format :fasta
                                 :callback #'make-seq-from-alist)))
      (is (eql 'iupac-dna-sequence (type-of s)))
      (is (eql *iupac-dna* (alphabet-of s)))
      (is (= 210 (length-of s)))
      (is (string= "Test1" (name-of s))))))

(test read-bio-sequence/multifasta/dna
  (with-open-file (stream (merge-pathnames "data/simple-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((line-buffer (make-line-buffer stream)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence line-buffer
                                    :alphabet :dna
                                    :format :fasta
                                    :callback #'make-seq-from-alist)))
          (is (eql 'simple-dna-sequence (type-of s)))
          (is (eql *dna* (alphabet-of s)))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (name-of s)))))
      (is (null (read-bio-sequence line-buffer
                                   :alphabet :dna
                                   :format :fasta)))))
  (with-open-file (stream (merge-pathnames "data/iupac-dna2.fa")
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let* ((line-buffer (make-line-buffer stream)))
      (dotimes (n 2)
        (let ((s (read-bio-sequence line-buffer
                                    :alphabet :dna
                                    :ambiguity :iupac
                                    :format :fasta
                                    :callback #'make-seq-from-alist)))
          (is (eql 'iupac-dna-sequence (type-of s)))
          (is (eql *iupac-dna* (alphabet-of s)))
          (is (= 280 (length-of s)))
          (is (string= (format nil "Test~a" (1+ n)) (name-of s)))))
      (is (null (read-bio-sequence line-buffer
                                   :alphabet :dna
                                   :format :fasta))))))
