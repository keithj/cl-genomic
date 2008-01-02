
(in-package :bio-sequence)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defclass alphabet ()
  ((name :initform nil
         :initarg :name
         :reader name-of
         :documentation "The alphabet name.")
   (tokens :initform ""
           :initarg :tokens
           :reader tokens-of
           :documentation "The set of member tokens of the
alphabet."))
    (:documentation "Alphabets are sets of tokens."))

  (defvar *dna*
    (make-instance 'alphabet
                   :name 'dna
                   :tokens (make-array 4
                                       :element-type 'base-char
                                       :initial-contents "acgt")))
  (defvar *rna*
    (make-instance 'alphabet
                   :name 'rna
                   :tokens (make-array 4
                                       :element-type 'base-char
                                       :initial-contents "acgu")))
  (defvar *iupac-dna*
    (make-instance 'alphabet
                   :name 'iupac-dna
                   :tokens
                   (make-array 15
                               :element-type 'base-char
                               :initial-contents "acgtrykmswbdhvn")))
  (defvar *iupac-rna*
    (make-instance 'alphabet
                   :name 'iupac-rna
                   :tokens
                   (make-array 15
                               :element-type 'base-char
                               :initial-contents "acgurykmswbdhvn"))))

(defclass bio-sequence ()
  ((alphabet :initarg :alphabet
             :accessor alphabet-of
             :documentation "The alphabet whose tokens comprise the
sequence.")
   (token-seq :initform nil
              :initarg :token-seq
              :accessor token-seq-of
              :documentation "The residue tokens of the sequence."))
  (:documentation "A biological sequence comprising an ordered string
of tokens from a specified alphabet."))

(defclass encoded-mixin ()
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
  ((metric :initform nil
           :initarg :metric
           :reader metric-of
           :documentation "A description of the quality metric
measured by the quality values. For example, p-value, Phred score or
Illumina score. This should be changed to a controlled vocabulary or
enumeration.")
   (quality :initform nil
            :initarg :quality
            :accessor quality-of
            :documentation "The array of quality values which should
be the same length as the array of residue tokens."))
  (:documentation "A mixin with support for bio-sequences that have a
numeric quality value for each residue."))

(defclass nucleic-acid-sequence (bio-sequence)
  ()
  (:documentation "A nucleic acid sequence comprising an ordered
string of tokens from a specified alphabet."))

(defclass dna-sequence (nucleic-acid-sequence)
  ()
  (:documentation "A DNA sequence."))

(defclass rna-sequence (nucleic-acid-sequence)
  ()
  (:documentation "An RNA sequence."))

(defclass simple-dna-sequence (dna-sequence encoded-mixin)
  ((alphabet :initform *dna*
             :allocation :class)
   (encoder :initform #'encode-dna-2bit
            :allocation :class)
   (decoder :initform #'decode-dna-2bit
            :allocation :class))
  (:documentation "A DNA sequence comprising unambiguous bases T, C, A
and G."))

(defclass iupac-dna-sequence (dna-sequence encoded-mixin)
  ((alphabet :initform *iupac-dna*
             :allocation :class)
   (encoder :initform #'encode-dna-4bit
            :allocation :class)
   (decoder :initform #'decode-dna-4bit
            :allocation :class))
  (:documentation "A DNA sequence comprising IUPAC ambiguity bases."))

(defclass simple-rna-sequence (rna-sequence encoded-mixin)
  ((alphabet :initform *rna*
             :allocation :class)
   (encoder :initform #'encode-rna-2bit
            :allocation :class)
   (decoder :initform #'decode-rna-2bit
            :allocation :class))
  (:documentation "An RNA sequence comprising unambiguous bases U, C,
A and G."))

(defclass iupac-rna-sequence (rna-sequence encoded-mixin)
  ((alphabet :initform *iupac-rna*
             :allocation :class)
   (encoder :initform #'encode-rna-4bit
            :allocation :class)
   (decoder :initform #'decode-rna-4bit
            :allocation :class))
  (:documentation "An RNA sequence comprising IUPAC ambiguity
bases."))

(defclass simple-dna-quality-sequence (simple-dna-sequence quality-mixin)
  ())

(defclass iupac-dna-quality-sequence (iupac-dna-sequence quality-mixin)
  ())
