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

(in-package :bio-sequence)

;;; identity-mixin generics

(defgeneric anonymousp (object)
  (:documentation "Returns T if OBJECT is anonymous, that is, has an
identifier of NIL."))


;;; bio-sequence generics

(defgeneric simplep (bio-sequence alphabet-designator)
  (:documentation "Returns T if BIO-SEQUENCE contains no ambiguous
residues according to the alphabets specified by ALPHABET-DESIGNATOR,
or NIL otherwise."))

(defgeneric size-of (alphabet)
  (:documentation "Returns the number of tokens in ALPHABET as a
fixnum."))

(defgeneric memberp (alphabet char)
  (:documentation "Returns T if CHAR is a member token of ALPHABET, or
NIL otherwise."))

(defgeneric virtualp (bio-sequence)
  (:documentation "Returns T if BIO-SEQUENCE has no concrete token-seq
representation, merely a length, or NIL otherwise."))

(defgeneric length-of (bio-sequence)
  (:documentation "Returns the length of BIO-SEQUENCE."))

(defgeneric residue-of (bio-sequence index)
  (:documentation "Returns the residue at INDEX of BIO-SEQUENCE. As
these are interbase coordinates, a residue actually occupies a range
from INDEX to 1+ INDEX. However, as the end index of the range is
always 1+ the start index, it is implied for convenience. The first
base of a sequence is therefore index 0 to index 1, with 0 being the
argument passed to the RESIDUE-OF method."))

(defgeneric (setf residue-of) (value bio-sequence index)
  (:documentation "Sets the residue at INDEX of BIO-SEQUENCE to
VALUE."))

(defgeneric to-string (bio-sequence &key start end token-case)
  (:documentation "Returns the string representing BIO-SEQUENCE,
starting at the first residue, or index START, and continuing to the
last residue, or index END."))

(defgeneric reverse-sequence (bio-sequence)
  (:documentation "Returns a new bio-sequence with a token-seq that is
a reversed copy of that of BIO-SEQUENCE."))

(defgeneric nreverse-sequence (bio-sequence)
  (:documentation "Returns BIO-SEQUENCE, destructively modified such
that its token-seq is reversed."))

(defgeneric complement-sequence (nucleic-acid-sequence)
  (:documentation "Returns a new nucleic-acid-sequence with a
token-seq that is a complemented copy of that of BIO-SEQUENCE."))

(defgeneric ncomplement-sequence (nucleic-acid-sequence)
  (:documentation "Returns NUCLEIC-ACID-SEQUENCE, destructively
modified such that its token-seq is complemented."))

(defgeneric reverse-complement (nucleic-acid-sequence)
  (:documentation "Returns a new nucleic-acid-sequence that is the
reverse-complement of NUCLEIC-ACID-SEQUENCE."))

(defgeneric nreverse-complement (nucleic-acid-sequence)
  (:documentation "Returns NUCLEIC-ACID-SEQUENCE, destructively
modified such that its token-seq is reverse-complemented."))

(defgeneric subsequence (bio-sequence start &optional end)
  (:documentation "Subsequence creates a bio-sequence that is a copy
of the subsequence of BIO-SEQUENCE bounded by START and END."))

(defgeneric search-sequence (bio-sequence-1 bio-sequence-2
                             &key from-end start1 start2 end1 end2)
  (:documentation "Searches bio-sequence-2 for a subsequence that
matches bio-sequence-1."))

(defgeneric cumulative-lengths (ranges)
  (:documentation "Returns a list of integers which are the
cumulative total lengths of RANGES."))


;;; bio-sequence io generics

(defgeneric begin-object (parser)
  (:documentation "Signals to PARSER the beginning of a new
object. The parser may use this information to discard accumulated
state from the previous object."))

(defgeneric object-relation (parser relation value)
  (:documentation "Signals to PARSER that RELATION exits between the
current object and some VALUE."))

(defgeneric object-identity (parser identity)
  (:documentation "Signals to PARSER that the current object has
IDENTITY."))

(defgeneric object-description (parser description)
  (:documentation "Signals to PARSER that the current object has
DESCRIPTION."))

(defgeneric object-alphabet (parser alphabet)
  (:documentation "Signals to PARSER that the current object has
residues of ALPHABET."))

(defgeneric object-residues (parser residues)
  (:documentation "Signals to PARSER that the current object has
bio-sequence RESIDUES. Multiple calls to this method may be used to
accumulate sequence residues."))

(defgeneric object-quality (parser quality)
  (:documentation "Signals to parser that the current object has
bio-sequence residue QUALITY. Multiple calls to this method may be
used to accumulate quality information."))

(defgeneric end-object (parser)
  (:documentation "Signals to PARSER the end of the current
object. The method must return T if the accumulated state of the
current object is valid, or NIL otherwise. If PARSER is constructing a
CLOS object, the object must be returned by this method."))

(defgeneric make-input-fn (stream format &key alphabet parser
                           &allow-other-keys)
  (:documentation "Returns a function of zero arity that uses PARSER
to read a single bio-sequence of ALPHABET, in FORMAT, from STREAM. The
function must return T on success, ot NIL otherwise."))

(defgeneric make-output-fn (stream format &key token-case
                            &allow-other-keys)
  (:documentation "Returns a function of arity 1 that accepts a
bio-sequence and writes a representation in FORMAT to STREAM. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids."))

(defgeneric make-bio-sequence (parser)
  (:documentation "Returns a new bio-sequence created from the state
accumulated by PARSER."))

(defgeneric read-fasta-sequence (stream alphabet parser)
  (:documentation "Reads a single Fasta format record of ALPHABET from
STREAM using PARSER."))

(defgeneric write-fasta-sequence (bio-sequence stream &key token-case)
  (:documentation "Writes BIO-SEQUENCE to STREAM in Fasta format. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids."))

(defgeneric read-fastq-sequence (stream alphabet parser)
  (:documentation "Reads a single Fastq format record of ALPHABET from
STREAM using PARSER."))

(defgeneric write-fastq-sequence (bio-sequence stream &key token-case)
  (:documentation "Writes BIO-SEQUENCE to STREAM in Fastq format. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids."))
