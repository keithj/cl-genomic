
(in-package :bio-sequence)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defclass alphabet ()
    ((name :initarg :name
           :reader name-of
           :documentation "The alphabet name.")
     (tokens :initform ""
             :initarg :tokens
             :reader tokens-of
             :documentation "The set of member tokens of the
alphabet."))
    (:documentation "Alphabets are sets of tokens."))

  (defclass sequence-strand ()
    ((name :initarg :name
           :reader name-of
           :documentation "The strand name.")
     (token :initarg :token
            :reader token-of
          :documentation "The token representing the strand.")
     (number :initarg :number
             :reader number-of
             :documentation "The number representing the strand."))
    (:documentation "The strand of a nucleotide sequence."))

  (defvar *dna*
    (make-instance 'alphabet
                   :name :dna
                   :tokens (make-array 4
                                       :element-type 'base-char
                                       :initial-contents "acgt")))
  (defvar *rna*
    (make-instance 'alphabet
                   :name :rna
                   :tokens (make-array 4
                                       :element-type 'base-char
                                       :initial-contents "acgu")))
  (defvar *iupac-dna*
    (make-instance 'alphabet
                   :name :iupac-dna
                   :tokens
                   (make-array 15
                               :element-type 'base-char
                               :initial-contents "acgtrykmswbdhvn")))
  (defvar *iupac-rna*
    (make-instance 'alphabet
                   :name :iupac-rna
                   :tokens
                   (make-array 15
                               :element-type 'base-char
                               :initial-contents "acgurykmswbdhvn")))

  (defvar *forward-strand*
    (make-instance 'sequence-strand
                   :name :forward
                   :token #\+
                   :number 1))

  (defvar *reverse-strand*
    (make-instance 'sequence-strand
                   :name :reverse
                   :token #\-
                   :number -1))

  (defvar *without-strand*
    (make-instance 'sequence-strand
                   :name :unstranded
                   :token #\.
                   :number 0))
  (defvar *unknown-strand*
    (make-instance 'sequence-strand
                   :name :unknown
                   :token #\?
                   :number nil)))

(defclass sequenced-mixin ()
  ((token-seq :initform (error "a :token-seq argument is required")
              :initarg :token-seq
              :accessor token-seq-of
              :documentation "The residue tokens of the sequence."))
  (:documentation "A mixin with support for bio-sequences whose
residues have been determined directly by experiment."))

(defclass encoded-mixin (sequenced-mixin)
  ((encoder :initarg :encoder
            :reader encoder-of
            :documentation "A function that accepts a single argument
and returns its encoded sybmcol value.")
   (decoder :initarg :decoder
            :reader decoder-of
            :documentation "A function that accepts a single encoded
token datum argument and returns the decoded value."))
  (:documentation "A mixin with support for bio-sequences where the
residue tokens are encoded from characters to a more compact binary
format."))

(defclass quality-mixin ()
  ((metric :initform (error "a :metric argument is required")
           :initarg :metric
           :reader metric-of
           :documentation "A description of the quality metric
measured by the quality values. For example, p-value, Phred score or
Illumina score. This should be changed to a controlled vocabulary or
enumeration.")
   (quality :initform (error "a :quality argument is required")
            :initarg :quality
            :accessor quality-of
            :documentation "The array of quality values which should
be the same length as the array of residue tokens."))
  (:documentation "A mixin with support for bio-sequences that have a
numeric quality value for each residue."))

(defclass located-mixin ()
  ((locations :initform nil
              :initarg :locations
              :accessor locations-of
              :documentation "A list of sequence-location objects.")))

(defclass sequence-location ()
  ((sequence :initarg :sequence
             :reader sequence-of
             :documentation "The bio-sequence in which the location
exists.")
   (position :initarg :position
             :reader position-of
             :documentation "The logical position of the location in
the referenced bio-sequence. This position is with respect to the
referenced sequence and may be a negative or a positive value greater
than the length of the referenced sequence.")))

(defclass bio-sequence (ontology-instance-mixin located-mixin)
  ((alphabet :initarg :alphabet
             :accessor alphabet-of
             :documentation "The alphabet whose tokens comprise the
sequence.")
   (length :initarg :length
           :accessor length-of
           :documentation "The logical length of the sequence in
residues. This value is not required to be supported by concrete
sequence data."))
  (:documentation "A biological sequence comprising tokens from a
specified alphabet. Its position relative to other bio-sequences may
be given by specifying the sequence-locations of its start in those
sequences."))

(defclass nucleic-acid-sequence (bio-sequence)
  ()
  (:documentation "A logical nucleic acid sequence."))

(defclass dna-sequence (nucleic-acid-sequence)
  ()
  (:documentation "A logical DNA sequence."))

(defclass rna-sequence (nucleic-acid-sequence)
  ()
  (:documentation "A logical RNA sequence."))

(defclass iupac-dna-sequence (encoded-mixin dna-sequence)
  ((alphabet :initform *iupac-dna*
             :allocation :class)
   (encoder :initform #'encode-dna-4bit
            :allocation :class)
   (decoder :initform #'decode-dna-4bit
            :allocation :class))
  (:documentation "A concrete DNA sequence comprising IUPAC ambiguity
bases."))

(defclass simple-dna-sequence (encoded-mixin dna-sequence)
  ((alphabet :initform *dna*
             :allocation :class)
   (encoder :initform #'encode-dna-2bit
            :allocation :class)
   (decoder :initform #'decode-dna-2bit
            :allocation :class))
  (:documentation "A concrete DNA sequence comprising unambiguous
bases T, C, A and G."))

(defclass iupac-rna-sequence (encoded-mixin rna-sequence)
  ((alphabet :initform *iupac-rna*
             :allocation :class)
   (encoder :initform #'encode-rna-4bit
            :allocation :class)
   (decoder :initform #'decode-rna-4bit
            :allocation :class))
  (:documentation "A concrete RNA sequence comprising IUPAC ambiguity
bases."))

(defclass simple-rna-sequence (encoded-mixin rna-sequence)
  ((alphabet :initform *rna*
             :allocation :class)
   (encoder :initform #'encode-rna-2bit
            :allocation :class)
   (decoder :initform #'decode-rna-2bit
            :allocation :class))
  (:documentation "A concrete RNA sequence comprising unambiguous
bases U, C, A and G."))

(defclass simple-dna-quality-sequence (simple-dna-sequence quality-mixin)
  ())

(defclass iupac-dna-quality-sequence (iupac-dna-sequence quality-mixin)
  ())
