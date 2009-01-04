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
(defgeneric simplep (bio-sequence)
  (:documentation "Returns T if BIO-SEQUENCE contains no ambiguous
residues, or NIL otherwise."))

(defgeneric ambiguousp (bio-sequence)
  (:documentation "Returns T if BIO-SEQUENCE contains ambiguous
residues, or NIL otherwise."))

(defgeneric size-of (alphabet)
  (:documentation "Returns the number of tokens in ALPHABET as a
fixnum."))

(defgeneric memberp (char alphabet)
  (:documentation "Returns T if CHAR is a member token of ALPHABET, or
NIL otherwise."))

(defgeneric subsumesp (token1 token2 alphabet)
  (:documentation ""))

(defgeneric token-index (token alphabet)
  (:documentation ""))

(defgeneric enum-ambiguity (char alphabet)
  (:documentation "Returns a list of the ambiguity characters
represented by CHAR in ALPHABET."))

(defgeneric strand-designator-p (object)
  (:documentation ""))

(defgeneric decode-strand (strand-designator &key strict)
  (:documentation ""))

(defgeneric invert-strand (sequence-strand)
  (:documentation ""))

(defgeneric match-strand (sequence-strand sequence-strand)
  (:documentation ""))

(defgeneric forward-strand-p (sequence-strand)
  (:documentation ""))

(defgeneric reverse-strand-p (sequence-strand)
  (:documentation ""))

(defgeneric unknown-strand-p (sequence-strand)
  (:documentation ""))

(defgeneric element-of (bio-sequence index)
  (:documentation "Returns the element at INDEX of BIO-SEQUENCE.  As
these are interbase coordinates, a residue actually occupies a range
from INDEX to 1+ INDEX. However, as the end index of the range is
always 1+ the start index, it is implied for convenience. The first
base of a sequence is therefore index 0 to index 1, with 0 being the
argument passed to the INDEX-OF method."))

(defgeneric (setf element-of) (value bio-sequence index)
  (:documentation "Sets the residue at INDEX of BIO-SEQUENCE to
VALUE."))

(defgeneric virtualp (bio-sequence)
  (:documentation "Returns T if BIO-SEQUENCE has no concrete token
sequence representation, merely a length, or NIL otherwise."))

(defgeneric single-stranded-p (nucleic-acid-sequence)
  (:documentation ""))

(defgeneric double-stranded-p (nucleic-acid-sequence)
  (:documentation ""))

(defgeneric num-gaps-of (bio-sequence &key start end)
  (:documentation ""))

(defgeneric residue-frequencies (bio-sequence alphabet)
  (:documentation ""))

(defgeneric coerce-sequence (bio-sequence result-type &key start end)
  (:documentation "Coerces all or part of BIO-SEQUENCE to
RESULT-TYPE. RESULT-TYPE may be 'string or another bio-sequence
class. For example an object of class {defclass dna-sequence} may be
coerced to an {defclass rna-sequence} and vice versa. The coerced
sequence starts at the first residue, or index START, and continuing
to the last residue, or index END."))

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

(defgeneric translate (nucleic-acid-sequence genetic-code
                       &key start end initiator-codon partial-codon))

(defgeneric translate-codon (codon genetic-code &key initiator)
  (:documentation "Returns an amino acid for CODON using
  GENETIC-CODE."))

(defgeneric start-codon-p (codon genetic-code)
  (:documentation "Returns T if CODON may be a start codon in
  GENETIC-CODE."))

(defgeneric term-codon-p (codon genetic-code)
  (:documentation "Returns T if CODON may be a terminator in
  GENETIC-CODE."))

(defgeneric residue-position (character bio-sequence
                              &key from-end test test-not start end)
  (:documentation ""))

(defgeneric quality-position (quality bio-sequence
                              &key from-end test test-not start end)
  (:documentation ""))


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

(defgeneric make-seq-input (stream format &key alphabet parser
                            &allow-other-keys)
  (:documentation "Returns a generator function of arity 1 that uses
PARSER to read a single bio-sequence of ALPHABET, in FORMAT, from
STREAM. The standard generator interface functions, CURRENT, NEXT and
HAS-MORE-P may be used in operations on the returned generator."))

(defgeneric make-seq-output (stream format &key token-case
                             &allow-other-keys)
  (:documentation "Returns a consumer function of arity 1 that accepts
a bio-sequence and writes a representation in FORMAT to STREAM. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids. The standard consumer interface function
CONSUME may be used in operations on the returned consumer."))

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

(defgeneric split-sequence-file (filespec format pathname-gen
                                 &key chunk-size)
  (:documentation "Splits sequence file identified by FILESPEC into
automatically named files, each containing up to CHUNK-SIZE
records. The new file names are created by the function
PATHNAME-GEN. See iou:pathname-generator and iou:pathname-extender
."))

(defgeneric convert-sequence-file (in-filespec in-format
                                   out-filespec out-format)
  (:documentation "Converts the sequence data in the file identified
by IN-FILESPEC in format IN-FORMAT, to OUT-FORMAT, writing the data to
a new file identified by OUT-FILESPEC. Returns the number of records
converted."))

(defgeneric align-local (seqm seqn subst-fn &key gap-open gap-extend
                         band-centre band-width alignment)
  (:documentation "Performs Smith Waterman sequence alignment using
affine gap penalties. The affine gap scoring is expressed as
penalties, that is to say GAP-OPEN and GAP-EXTEND should be negative
values. The alignment may be banded to prune the search space.

Arguments:

- seqm (object): A sequence to be aligned.
- seqn (object): A sequence to be aligned.
- subst-fn (function): A substitution function that accepts two
sequence elements as arguments and returns a single-float substitution
score.

Key:

- gap-open (single-float): The gap opening score, a negative value.

- gap-extend (single-float): The gap extension score, a negative
  value.

- band-centre (fixnum): The band centre for banded alignments. This
defaults to 0, the main diagonal. The desired band may be calculated
by subtracting a j coordinate from its corresponding i coordinate.
- band-width (fixnum): The band width for banded alignments. This
defaults to most-positive-fixnum so that the search space is not
pruned.

- alignment (generalized boolean): T if an alignment is to be
  calculated.

Returns:

- A single-float alignment score.
- An alignment object, or NIL."))

(defgeneric align-local-ksh (seqm seqn subst-fn &key gap-open gap-extend k
                             alignment)
  (:documentation "Performs Smith Waterman sequence alignment using
affine gap penalties with a kmer seeding heuristic. The affine gap
scoring is expressed as penalties, that is to say GAP-OPEN and
GAP-EXTEND should be negative values. The alignment is seeded with
kmers of length K and the alignment is banded to include all such kmers.

Arguments:

- seqm (object): A sequence to be aligned.
- seqn (object): A sequence to be aligned.
- subst-fn (function): A substitution function that accepts two
sequence elements as arguments and returns a single-float substitution
score.

Key:

- gap-open (single-float): The gap opening penalty.
- gap-extend (single-float): The gap extension penalty.
- k (fixnum): The seed kmer length.

- alignment (generalized boolean): T if an alignment is to be
  calculated.

Returns:

- A single-float alignment score.
- An alignment object, or NIL."))


(defgeneric  begin-term (parser)
  (:documentation ""))

(defgeneric  begin-typedef (parser)
  (:documentation ""))

(defgeneric  begin-instance (parser)
  (:documentation ""))

(defgeneric  end-section (parser)
  (:documentation ""))

(defgeneric  tag-value (parser tag value)
  (:documentation ""))

