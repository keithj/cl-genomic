
(in-package :bio-sequence)

(defclass alphabet ()
  ((name :initform nil
         :initarg :name
         :reader name-of
         :documentation "The alphabet name.")
   (symbols :initform ""
            :initarg :symbols
            :reader symbols-of
            :documentation "The set of member symbols of the
alphabet."))
  (:documentation "Alphabets are sets of symbols."))

(defvar *dna*
  (make-instance 'alphabet
                 :name 'dna
                 :symbols (make-array 4
                                      :element-type 'base-string
                                      :initial-contents "acgt")))

(defvar *rna*
  (make-instance 'alphabet
                 :name 'rna
                 :symbols (make-array 4
                                      :element-type 'base-string
                                      :initial-contents "acgu")))

(defvar *iupac-dna*
   (make-instance 'alphabet
                 :name 'iupac-dna
                 :symbols
                 (make-array 15 :element-type 'base-string
                             :initial-contents "acgtrykmswbdhvn")))

(defvar *iupac-rna*
   (make-instance 'alphabet
                 :name 'iupac-rna
                 :symbols
                 (make-array 15 :element-type 'base-string
                             :initial-contents "acgurykmswbdhvn")))

(defclass bio-sequence ()
  ((alphabet :initarg :alphabet
             :accessor alphabet-of
             :documentation "The alphabet whose symbols comprise the
sequence.")
   (symbol-seq :initform nil
               :initarg :symbol-seq
               :accessor symbol-seq-of
               :documentation "The residue symbols of the
sequence."))
  (:documentation "A biological sequence comprising an ordered string
of symbols from a specified alphabet."))

(defclass encoded-mixin ()
  ((encoder :initarg :encoder
            :reader encoder-of
            :documentation "A function that accepts a single argument
and returns its encoded sysmcol value.")
   (decoder :initarg :decoder
            :reader decoder-of
            :documentation "A function that accepts a single encoded
symbol datum argument and returns the decoded value."))
  (:documentation "The encoding used to convert an alphabet's
symbols."))

(defclass nucleic-acid-sequence (bio-sequence)
  ()
  (:documentation "A nucleic acid sequence comprising an ordered
string of symbols from a specified alphabet."))

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
  (:documentation "An RNA sequence comprising unambiguous bases U, C, A
and G."))

(defclass iupac-rna-sequence (rna-sequence encoded-mixin)
  ((alphabet :initform *iupac-rna*
             :allocation :class)
   (encoder :initform #'encode-rna-4bit
            :allocation :class)
   (decoder :initform #'decode-rna-4bit
            :allocation :class))
  (:documentation "An RNA sequence comprising IUPAC ambiguity
bases."))
