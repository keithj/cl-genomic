
(defpackage bio-sequence
  (:use #:common-lisp #:cl-io-utilities #:cl-gp-utilities)
  (:export
   ;; Specials
   #:*dna*
   #:*rna*
   #:*iupac-dna*
   #:*iupac-rna*
   ;; Functions
   #:make-dna-seq
   #:make-rna-seq
   #:make-dna-quality-seq
   #:phred-quality
   #:encode-phred-quality
   #:decode-phred-quality
   #:illumina-quality
   #:encode-illumina-quality
   #:decode-illumina-quality
   #:illumina-to-phred-quality
   ;; Classes
   #:alphabet
   #:bio-sequence
   #:nucleic-acid-sequence
   #:dna-sequence
   #:rna-sequence
   #:simple-dna-sequence
   #:simple-rna-sequence
   #:iupac-dna-sequence
   #:iupac-rna-sequence
   #:simple-dna-quality-sequence
   #:iupac-dna-quality-sequence
   ;; Generic functions
   #:name-of
   #:tokens-of
   #:alphabet-of
   #:residue-of
   #:length-of
   #:copy-sequence
   #:to-string
   #:metric-of
   #:quality-of))
