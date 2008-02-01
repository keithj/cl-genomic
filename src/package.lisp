
(defpackage bio-sequence
  (:use #:common-lisp #:cl-io-utilities #:cl-gp-utilities
        #:trivial-gray-streams #:split-sequence)
  (:export
   ;; Specials
   #:*dna*
   #:*rna*
   #:*iupac-dna*
   #:*iupac-rna*

   #:*default-knowledgebase*
   ;; Conditions
   #:knowledgebase-error

   ;; Functions
   #:read-bio-sequence
   #:make-seq
   #:make-quality-seq
   #:make-seq-from-alist
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

   #:knowledgebase
   #:frame
   #:slot
   #:reflexive-mixin
   #:transitive-mixin
   #:inverse-mixin
   #:part-of
   #:has-part
   #:subsequence-of
   #:has-subsequence
   #:instance

   ;; Generic functions
   #:identity-of
   #:name-of
   #:tokens-of
   #:alphabet-of
   #:residue-of
   #:length-of
   #:copy-sequence
   #:to-string
   #:metric-of
   #:quality-of
   #:read-bio-sequence

   #:contains-frame-p
   #:find-frame
   #:add-frame
   #:remove-frame
   #:contains-slot-p
   #:find-slot
   #:slots-of
   #:add-slot
   #:remove-slot
   #:domain-of
   #:range-of
   #:value-of
   #:slot-value-of
   ))
