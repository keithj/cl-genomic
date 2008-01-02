
(in-package :cl-bio-system)

(fiveam:def-suite testsuite
    :description "The test suite.")


(in-package :cl-bio-test)

(fiveam:in-suite cl-bio-system:testsuite)

;; Test alphabet
(test name-of/alphabet
  (let ((name "test name"))
    (is (string= name
                 (name-of (make-instance 'alphabet
                                         :name name))))))

(test tokens-of/alphabet
  (let ((tokens "abcd"))
    (is (string= tokens
                 (tokens-of (make-instance 'alphabet
                                           :tokens tokens))))))

;; Test encoding/decoding
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

;; Test constructors
(test make-bio-seq
  (let ((seqs (list (make-dna-seq "tagc" :ambiguity nil)
                    (make-dna-seq "nnnn" :ambiguity :iupac)
                    (make-rna-seq "uagc" :ambiguity nil)
                    (make-rna-seq "nnnn" :ambiguity :iupac)))
        (class-names (list 'simple-dna-sequence
                           'iupac-dna-sequence
                           'simple-rna-sequence
                           'iupac-rna-sequence)))
    (mapcar #'(lambda (seq class-name)
                (is (eql class-name (class-name (class-of seq)))))
            seqs class-names)))

(test make-quality-seq
  (let ((seqs (list (make-dna-quality-seq
                     "agaatattctgaccccagttactttcaaga"
                     "<<<<<<<<<<<<<<<<<<<<<735513;3<"
                     :ambiguity nil :metric :phred)
                    (make-dna-quality-seq
                     "ntgccaaaaaatagaaaagtcancgatatt"
                     "<<<<<<<<<<<<<<<<<8<<<<<<5<3.5:"
                     :ambiguity :iupac :metric :phred)))
        (class-names (list 'simple-dna-quality-sequence
                           'iupac-dna-quality-sequence)))
    (mapcar #'(lambda (seq class-name)
                (is (eql class-name (class-name (class-of seq)))))
            seqs class-names)))

;; Test bio-sequence methods
(test length-of/bio-sequence
  (let ((len 10))
    (is (= len (length-of (make-dna-seq "aaaaaaaaaa"
                                        :ambiguity nil))))
    (is (= len (length-of (make-rna-seq "aaaaaaaaaa"
                                        :ambiguity nil))))))

(test residue-of/dna-sequence
  (let ((residues "tttt")
        (seq (make-dna-seq "aaaa" :ambiguity nil)))
    (dotimes (n (length residues))
      (setf (residue-of seq n) (aref residues n))
      (is (char= (residue-of seq n) (aref residues n))))))

(test residue-of/rna-sequence
  (let ((residues "uuuu")
        (rna-seq (make-rna-seq "aaaa" :ambiguity nil)))
    (dotimes (n (length residues))
      (setf (residue-of rna-seq n) (aref residues n))
      (is (char= (residue-of rna-seq n) (aref residues n))))))

(test copy-sequence/bio-sequence
  (let ((seqs (list (make-dna-seq "tagc" :ambiguity nil)
                    (make-dna-seq "nnnn" :ambiguity :iupac)
                    (make-rna-seq "uagc" :ambiguity nil)
                    (make-rna-seq "nnnn" :ambiguity :iupac))))
    (mapcar #'(lambda (seq)
                (let ((copy (copy-sequence seq)))
                  (is (= (length-of copy) (length-of seq)))
                  (dotimes (n 4)
                    (is (char= (residue-of copy n)
                               (residue-of seq n))))))
            seqs)))

(test to-string/dna-sequence
  (let* ((residues "acgt")
         (seq (make-dna-seq residues :ambiguity nil)))
    ;; no args
    (is (string= residues (to-string seq)))
    ;; optional arg start
    (dotimes (n 4)
      (is (string= (subseq residues n) (to-string seq n))))
    (dotimes (n 4)
      (is (string= (subseq residues 0 n) (to-string seq 0 n))))))
