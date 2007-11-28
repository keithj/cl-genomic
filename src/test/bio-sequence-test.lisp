
(in-package :cl-bio-test)

;; Test alphabet
(define-test name-of/alphabet
  (let ((name "test name"))
    (assert-equal name
                  (name-of (make-instance 'alphabet
                                          :name name)))))

(define-test symbols-of/alphabet
  (let ((symbols "abcd"))
    (assert-equal symbols
                  (symbols-of (make-instance 'alphabet
                                             :symbols symbols)))))

;; Test encoding/decoding
(define-test encode/decode-dna-2bit
  (let ((symbols "tagc"))
    (loop for sym across symbols
         do (assert-equal
             sym (bio-sequence::decode-dna-2bit
                  (bio-sequence::encode-dna-2bit sym))))))

(define-test encode/decode-rna-2bit
  (let ((symbols "uagc"))
    (loop for sym across symbols
         do (assert-equal
             sym (bio-sequence::decode-rna-2bit
                  (bio-sequence::encode-rna-2bit sym))))))

(define-test encode/decode-dna-4bit
  (let ((symbols "tagcrykmswbdhvn"))
    (loop for sym across symbols
       do (assert-equal
           sym (bio-sequence::decode-dna-4bit
                (bio-sequence::encode-dna-4bit sym))))))

(define-test encode/decode-rna-4bit
  (let ((symbols "uagcrykmswbdhvn"))
    (loop for sym across symbols
         do (assert-equal
             sym (bio-sequence::decode-rna-4bit
                  (bio-sequence::encode-rna-4bit sym))))))

;; Test constructors
(define-test make-bio-seq
  (let ((seqs (list (make-dna-seq "tagc" :ambiguity nil)
                    (make-dna-seq "nnnn" :ambiguity :iupac)
                    (make-rna-seq "uagc" :ambiguity nil)
                    (make-rna-seq "nnnn" :ambiguity :iupac)))
        (class-names (list 'simple-dna-sequence
                           'iupac-dna-sequence
                           'simple-rna-sequence
                           'iupac-rna-sequence)))
    (mapcar #'(lambda (seq class-name)
                (assert-eql class-name (class-name (class-of seq))))
            seqs class-names)))

;; Test bio-sequence methods
(define-test length-of/bio-sequence
  (let ((len 10))
    (assert-eq len (length-of (make-dna-seq "aaaaaaaaaa"
                                            :ambiguity nil)))
    (assert-eq len (length-of (make-rna-seq "aaaaaaaaaa"
                                            :ambiguity nil)))))

(define-test residue-of/dna-sequence
  (let ((residues "tttt")
        (seq (make-dna-seq "aaaa" :ambiguity nil)))
    (dotimes (n (length residues))
      (setf (residue-of seq n) (aref residues n))
      (assert-eq (residue-of seq n) (aref residues n)))))

(define-test residue-of/rna-sequence
  (let ((residues "uuuu")
        (rna-seq (make-rna-seq "aaaa" :ambiguity nil)))
    (dotimes (n (length residues))
      (setf (residue-of rna-seq n) (aref residues n))
      (assert-eq (residue-of rna-seq n) (aref residues n)))))

(define-test copy-sequence/bio-sequence
  (let ((seqs (list (make-dna-seq "tagc" :ambiguity nil)
                    (make-dna-seq "nnnn" :ambiguity :iupac)
                    (make-rna-seq "uagc" :ambiguity nil)
                    (make-rna-seq "nnnn" :ambiguity :iupac))))
    (mapcar #'(lambda (seq)
                (let ((copy (copy-sequence seq)))
                  (assert-eq (length-of copy) (length-of seq))
                  (dotimes (n 4)
                    (assert-eq (residue-of copy n)
                               (residue-of seq n)))))
            seqs)))

(define-test to-string/dna-sequence
  (let* ((residues "acgt")
         (seq (make-dna-seq residues :ambiguity nil)))
    ;; no args
    (assert-equal residues (to-string seq))
    ;; optional arg start
    (dotimes (n 4)
      (assert-equal (subseq residues n) (to-string seq n)))
    (dotimes (n 4)
      (assert-equal (subseq residues 0 n) (to-string seq 0 n)))))
