;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
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
  (:documentation "Returns T if TOKEN1 subsumes TOKEN2 using IUPAC
ambiguities, where both tokens belong to ALPHABET."))

(defgeneric token-index (token alphabet)
  (:documentation "Returns the integer position of TOKEN in
ALPHABET."))

(defgeneric  random-token-of (alphabet)
  (:documentation "Returns a random token from ALPHABET, with an
equal probability of selecting each token."))

(defgeneric enum-ambiguity (char alphabet)
  (:documentation "Returns a list of the ambiguity characters
represented by CHAR in ALPHABET."))

(defgeneric strand-designator-p (object)
  (:documentation "Returns T if OBJECT is one of the supported nucleic
acid sequence strand designators, or NIL otherwise.

The strand designators are:

                   String  Character     Symbol    Integer
Forward strand          +          +   :forward          1
Reverse strand          -          -   :reverse         -1
Unknown strand          ?          ?   :unknown          0"))

(defgeneric decode-strand (strand-designator &key strict)
  (:documentation "Returns a {defclass sequence-strand} given
STRAND-DESIGNATOR, or NIL if the designator is not recognised. If
STRICT is set to T an invalid designator will raise an error."))

(defgeneric strand= (sequence-strand1 sequence-strand2)
  (:documentation "Returns T if SEQUENCE-STRAND1 can be shown to be
the same as SEQUENCE-STRAND2, or NIL otherwise. If either strand is
unknown, the result will be NIL."))

(defgeneric strand-equal (sequence-strand1 sequence-strand2)
  (:documentation "Returns T if SEQUENCE-STRAND1 may be the same as
SEQUENCE-STRAND2, or NIL otherwise. If either strand is unknown, the
result will be T."))

(defgeneric complement-strand (sequence-strand)
  (:documentation "Returns the complement strand to
SEQUENCE-STRAND. The complement of the unknown strand is the unknown
strand."))

(defgeneric complementp (sequence-strand1 sequence-strand2)
  (:documentation "Returns T if SEQUENCE-STRAND1 is the complement of
SEQUENCE-STRAND2, or NIL otherwise."))

(defgeneric forward-strand-p (sequence-strand)
  (:documentation "Returns T if SEQUENCE-STRAND is the forward strand,
or NIL otherwise."))

(defgeneric reverse-strand-p (sequence-strand)
  (:documentation "Returns T if SEQUENCE-STRAND is the reverse strand,
or NIL otherwise."))

(defgeneric unknown-strand-p (sequence-strand)
  (:documentation "Returns T if SEQUENCE-STRAND is the unknown strand,
or NIL otherwise."))

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
  (:documentation "Returns T if NUCLEIC-ACID-SEQUENCE is
single-stranded, or NIL otherwise."))

(defgeneric double-stranded-p (nucleic-acid-sequence)
  (:documentation "Returns T if NUCLEIC-ACID-SEQUENCE is
double-stranded, or NIL otherwise."))

(defgeneric num-gaps-of (bio-sequence &key start end)
  (:documentation ""))

(defgeneric residue-frequencies (bio-sequence)
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
reverse-complement of NUCLEIC-ACID-SEQUENCE. If the argument is a
{defclass na-sequence-interval} this is accomplished by returning a
copy at the corresponding position on the complementary strand of the
reference sequence."))

(defgeneric nreverse-complement (nucleic-acid-sequence)
  (:documentation "Returns NUCLEIC-ACID-SEQUENCE, destructively
modified such that its token-seq is reverse-complemented. If the
argument is a {defclass na-sequence-interval} this is accomplished by
moving it to the corresponding position on the complementary strand of
the reference sequence."))

(defgeneric subsequence (bio-sequence start &optional end)
  (:documentation "Subsequence creates a bio-sequence that is a copy
of the subsequence of BIO-SEQUENCE bounded by START and END."))

(defgeneric search-sequence (bio-sequence1 bio-sequence2
                             &key from-end start1 start2 end1 end2)
  (:documentation "Searches BIO-SEQUENCE2 for a subsequence that
matches BIO-SEQUENCE1."))

(defgeneric hamming-distance (bio-sequence1 bio-sequence2
                              &key start1 start2 end1 end2)
  (:documentation "Returns the Hamming distance between BIO-SEQUENCE1
and BIO-SEQUENCE2. Sequences or sequence ranges must be the same
length."))

(defgeneric hamming-search (bio-sequence1 bio-sequence2
                            &key start1 start2 end1 end2 max-distance)
  (:documentation "Searches BIO-SEQUENCE2 for a subsequence that
is up to MAX-DISTANCE Hamming distance from BIO-SEQUENCE1."))

(defgeneric cumulative-lengths (ranges)
  (:documentation "Returns a list of integers which are the
cumulative total lengths of RANGES."))

(defgeneric translate (nucleic-acid-sequence genetic-code
                       &key start end initiator-codon partial-codon)
  (:documentation "Translates a region of NUCLEIC-ACID-SEQUENCE,
bounded by START and END, using GENETIC-CODE. If INITIATOR-CODON is T,
the region is must begin with a recognised initiator codon, according
to GENETIC-CODE. If PARTIAL-CODON is T, the last codon in the region
may be partial and is padded with fully ambiguous bases at its end."))

(defgeneric translate-codon (codon genetic-code &key initiator)
  (:documentation "Returns an amino acid for CODON using
GENETIC-CODE. An implementation should consider ambiguities in CODON,
effectively enumerating all possible unambiguous codons from it,
translating them and unifying the result into a single, possibly
ambiguous amino acid."))

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

(defgeneric seguid (bio-sequence &key start end)
  (:documentation "Returns the SEquence Globally Unique IDentifier or
SEGUID string for BIO-SEQUENCE, between START and END. A SEGUID is
essentially a base64-encoded SHA-1 hash of the uppercase sequence
string. See PMID:16858731."))

(defgeneric md5 (bio-sequence &key start end)
  (:documentation "Returns the MD5 checksum of the uppercase sequence
string of BIO-SEQUENCE, between START and END."))

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

(defgeneric make-seq-input (source format &key alphabet parser
                            &allow-other-keys)
  (:documentation "Returns a generator function of arity 1 that uses
PARSER to read a single bio-sequence of ALPHABET, in FORMAT, from
SOURCE. The standard generator interface functions, CURRENT, NEXT and
HAS-MORE-P may be used in operations on the returned generator."))

(defgeneric make-seq-output (dest format &key token-case
                             &allow-other-keys)
  (:documentation "Returns a consumer function of arity 1 that accepts
a bio-sequence and writes a representation in FORMAT to DEST. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids. The standard consumer interface function
CONSUME may be used in operations on the returned consumer."))

(defgeneric make-bio-sequence (parser)
  (:documentation "Returns a new bio-sequence created from the state
accumulated by PARSER."))

(defgeneric make-interval-input (sequence)
  (:documentation "Returns a generator function that will return
elements of SEQUENCE, followed by a sentinel element."))

(defgeneric make-interval (bio-sequence &rest args)
  (:documentation "Returns a new interval on BIO-SEQUENCE."))

(defgeneric has-sequence-p (stream format &key alphabet)
  (:documentation "Returns T if a bio-sequence may be read from
STREAM."))

;;; Could all the read-XXXX-sequence and write-XXXX-sequence functions
;;; be replaced sensibly by read-bio-sequence and write-bio-sequence
;;; functions that dispatch on a format object?
(defgeneric read-pure-sequence (stream alphabet parser)
  (:documentation "Reads a single pure format record of ALPHABET from
STREAM using PARSER."))

(defgeneric write-pure-sequence (bio-sequence stream &key token-case)
  (:documentation "Writes BIO-SEQUENCE to STREAM in pure format. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids."))

(defgeneric read-raw-sequence (stream alphabet parser)
  (:documentation "Reads a single raw format record of ALPHABET from
STREAM using PARSER."))

(defgeneric write-raw-sequence (bio-sequence stream &key token-case)
  (:documentation "Writes BIO-SEQUENCE to STREAM in raw format. The
TOKEN-CASE keyword is used to override the default character case for
printing sequence residues, which is lowercase for DNA or RNA and
uppercase for amino acids."))

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
PATHNAME-GEN. See dxi:pathname-generator and dxi:pathname-extender."))

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
