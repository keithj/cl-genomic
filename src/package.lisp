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
  (:export
   ;; Specials
   #:*dna*
   #:*rna*

   #:*forward-strand*
   #:*reverse-strand*
   #:*without-strand*
   #:*unknown-strand*

   #:*default-knowledgebase*
   ;; Conditions
   #:knowledgebase-error

   ;; Functions
   #:complement-dna
   #:complement-rna
   #:find-alphabet
   #:read-bio-sequence
   #:write-bio-sequence
   #:make-quality-dna
   #:phred-quality
   #:encode-phred-quality
   #:decode-phred-quality
   #:illumina-quality
   #:encode-illumina-quality
   #:decode-illumina-quality
   #:illumina-to-phred-quality

   #:read-seq-datum
   #:write-seq-datum
   #:filter-seq-datum
   #:make-seq-datum
   #:make-quality-datum
   #:seq-datum-identity
   #:seq-datum-alphabet
   #:seq-datum-ambiguity
   #:seq-datum-token-seq
   #:seq-datum-length
   #:seq-datum-description
   #:seq-datum-quality
   #:make-seq-from-datum
   
   ;; Classes
   #:alphabet
   #:bio-sequence
   #:nucleic-acid-sequence
   #:dna-sequence
   #:rna-sequence
   #:dna-quality-sequence

   #:knowledgebase
   #:reflexive-mixin
   #:transitive-mixin
   #:inverse-mixin
   #:frame
   #:single-valued-slot
   #:set-valued-slot
   #:single-valued-inverse-slot
   #:set-valued-inverse-slot
   #:part-of
   #:has-part
   #:subsequence-of
   #:has-subsequence
   #:instance

   ;; Generic functions
   #:anonymousp
   #:identity-of
   #:name-of
   #:simplep
   #:virtualp
   #:size-of
   #:memberp
   #:tokens-of
   #:encoded-index-of
   #:decoded-index-of
   #:alphabet-of
   #:residue-of
   #:length-of
   #:subsequence
   #:reverse-sequence
   #:nreverse-sequence
   #:complement-sequence
   #:ncomplement-sequence
   #:reverse-complement
   #:nreverse-complement
   #:residue-frequencies
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
   #:slot-value-of))

(defpackage bio-sequence-user
  (:use #:common-lisp #:bio-sequence
        #:cl-io-utilities #:cl-gp-utilities
        #:trivial-gray-streams)
  (:nicknames #:bsu))
