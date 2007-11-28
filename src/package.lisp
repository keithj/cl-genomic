
(defpackage bio-sequence
  (:use #:common-lisp #:cl-gp-utilities)
  (:export
   ;; Specials
   #:*dna*
   #:*rna*
   #:*iupac-dna*
   #:*iupac-rna*
   ;; Functions
   #:make-dna-seq
   #:make-rna-seq
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
   ;; Generic functions
   #:name-of
   #:symbols-of
   #:alphabet-of
   #:residue-of
   #:length-of
   #:copy-sequence
   #:to-string
   #:from-string))
