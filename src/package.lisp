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

(defpackage bio-sequence
  (:use #:common-lisp #:cl-io-utilities #:cl-gp-utilities
        #:trivial-gray-streams #:split-sequence)
  (:nicknames #:bs)
  (:export
   ;; Types
   #:quality-score
   
   ;; Specials
   #:*sequence-print-limit*
   
   #:*forward-strand*
   #:*reverse-strand*
   #:*unknown-strand*

   #:*default-knowledgebase*
   
   ;; Conditions
   #:bio-sequence-error
   #:bio-sequence-warning
   #:bio-sequence-io-error
   #:bio-sequence-op-error

   ;; Functions
   #:make-dna
   #:make-rna
   #:make-dna-quality
   #:make-aa
   #:complement-dna
   #:complement-rna
   #:find-alphabet
   #:phred-quality
   #:encode-phred-quality
   #:decode-phred-quality
   #:illumina-quality
   #:encode-illumina-quality
   #:decode-illumina-quality
   #:illumina-to-phred-quality

   #:simple-dna-subst
   #:iupac-dna-subst
   #:blosum-50-subst
   
   #:split-fastq-file
   #:write-raw-fastq
   #:write-n-fastq

   ;; Classes
   #:alphabet
   #:sequence-strand
   #:bio-sequence
   #:na-sequence
   #:dna-sequence
   #:rna-sequence
   #:aa-sequence
   #:dna-quality-sequence
   #:interval
   #:bio-sequence-interval
   #:na-sequence-interval
   #:na-alignment-interval
   #:alignment

   #:bio-sequence-parser
   #:raw-sequence-parser
   #:simple-sequence-parser
   #:quality-sequence-parser
   #:virtual-sequence-parser

   ;; Generic functions
   #:anonymousp
   #:identity-of
   #:name-of
   #:symbol-of
   #:token-of
   #:number-of
   #:ambiguousp
   #:simplep
   #:virtualp
   #:strand-designator-p
   #:forward-strand-p
   #:reverse-strand-p
   #:unknown-strand-p
   #:decode-strand
   #:invert-strand
   #:match-strand
   #:num-strands-of
   #:single-stranded-p
   #:double-stranded-p
   #:num-gaps-of
   #:size-of
   #:memberp
   #:tokens-of
   #:explode-ambiguity
   #:encoded-index-of
   #:decoded-index-of
   #:alphabet-of
   #:element-of
   #:residue-of
   #:residues-of
   #:length-of
   #:metric-of
   #:quality-of
   #:upper-of
   #:lower-of
   #:strand-of
   #:reference-of
   #:subsequence
   #:reverse-sequence
   #:nreverse-sequence
   #:complement-sequence
   #:ncomplement-sequence
   #:reverse-complement
   #:nreverse-complement
   #:search-sequence
   #:residue-frequencies
   #:to-string
   #:invert-complement
   #:ninvert-complement

   #:aligned-of
   #:intervals-of
   #:aligned-length-of
   #:align-local
   #:align-local-ksh
   
   #:make-seq-input
   #:make-seq-output
   #:begin-object
   #:object-alphabet
   #:object-relation
   #:object-identity
   #:object-description
   #:object-residues
   #:object-quality
   #:end-object

   #:split-sequence-file
   #:convert-sequence-file
   
   #:parsed-alphabet
   #:parsed-identity
   #:parsed-description
   #:parsed-residues
   #:parsed-metric
   #:parsed-quality
   #:parsed-length
   #:parsed-raw)
  (:documentation "The BIO-SEQUENCE package provides basic support for
representing biological sequences. Concepts such as alphabets of
sequence residues (nucleic acid bases and amino acids), sequences of
residues (nucleic acids and polypeptides) and sequence strands are
represented as CLOS classes.

An event-based parsing interface enables reading of some common,
simple biological sequence file formats into primitive Lisp data
structures or CLOS instances."))

(defpackage bio-sequence-user
  (:use #:common-lisp #:bio-sequence
        #:cl-io-utilities #:cl-gp-utilities
        #:trivial-gray-streams)
  (:nicknames #:bsu))
