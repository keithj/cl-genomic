
(in-package :bio-sequence)

(deftype encoded-tokens (n)
  `(simple-array (unsigned-byte ,n) *))

(defvar *sequence-class-table*
  (pairlis '((:dna nil :default) (:rna nil :default)
             (:dna :iupac :default) (:rna :iupac :default)
             (:dna nil :quality)
             (:dna :iupac :quality))
           '(simple-dna-sequence simple-rna-sequence
             iupac-dna-sequence iupac-rna-sequence
             simple-dna-quality-sequence
             iupac-dna-quality-sequence)))

(defvar *sequence-print-limit* 50
  "Maximum length of sequence to be pretty-printed.")

(defmacro encode-token (value seq index token-seq-type)
  "Sets the residue token VALUE at INDEX in SEQ. The residue tokens
are encoded as an array of type TOKEN-SEQ-TYPE in SEQ."
  (let ((token-seq (gensym))
        (encoder (gensym)))
    `(let ((,token-seq (token-seq-of ,seq))
           (,encoder (encoder-of ,seq)))
       (declare (type ,token-seq-type ,token-seq)
                (type function ,encoder))
       (setf (aref ,token-seq ,index)
             (funcall ,encoder ,value)))))

(defmacro decode-token (seq index token-seq-type)
  "Returns the decoded residue token from INDEX in SEQ. The residue
tokens are encoded as an array of type TOKEN-SEQ-TYPE in SEQ."
  (let ((token-seq (gensym))
        (decoder (gensym)))
    `(let ((,token-seq (token-seq-of ,seq))
           (,decoder (decoder-of ,seq)))
       (declare (type ,token-seq-type ,token-seq)
                (type function ,decoder))
       (funcall ,decoder (aref ,token-seq ,index)))))

(defmacro decode-token-array (seq start end token-seq-type)
  "Returns a simple-base-string representing the token-seq of SEQ
from residues START to END, inclusive. The residues are encoded as
type TOKEN-SEQ-TYPE in the token-seq slot of SEQ."
  (let ((token-seq (gensym))
        (dest (gensym))
        (source-end (gensym))
        (dest-start (gensym))
        (decoder (gensym)))
    `(let ((,token-seq (token-seq-of ,seq))
           (,dest (make-string (- ,end ,start)
                               :element-type 'base-char))
           (,source-end (1- ,end))
           (,dest-start 0)
           (,decoder (decoder-of seq)))
       (declare (type ,token-seq-type ,token-seq)
                (type simple-base-string ,dest)
                (type array-index ,source-end ,dest-start)
                (type function ,decoder))
       (copy-array ,token-seq ,start ,source-end
                   ,dest ,dest-start ,decoder)
       ,dest)))

(defun encode-simple-seq (str encoder)
  "Encodes the tokens in simple-string STR with ENCODER and returns an
array of element type (unsigned-byte 2)."
  (declare (optimize (speed 3) (safety 1)))
  (declare (type simple-string str))
  (declare (type function encoder))
  (let ((token-seq (make-array (length str)
                                :element-type '(unsigned-byte 2))))
    (declare (type (encoded-tokens 2) token-seq))
    (copy-array str 0 (1- (length str))
                token-seq 0 encoder)
    token-seq))

(defun encode-iupac-seq (str encoder)
   "Encodes the tokens in simple-string STR with ENCODER and returns an
array of element type (unsigned-byte 4)."
  (declare (optimize (speed 3) (safety 1)))
  (declare (type simple-string str))
  (declare (type function encoder))
  (let ((token-seq (make-array (length str)
                                :element-type '(unsigned-byte 4))))
    (declare (type (encoded-tokens 4) token-seq))
    (copy-array str 0 (1- (length str))
                token-seq 0 encoder)
    token-seq))

(defun decode-quality (quality decoder)
  "Decodes the array QUALITY into a new array using function DECODER."
  (let ((quality-seq (make-array (length quality)
                                 :element-type '(unsigned-byte 8))))
    (copy-array quality 0 (1- (length quality))
                quality-seq 0 decoder)
    quality-seq))

(defun select-sequence-class (alphabet ambiguity quality)
  "Returns an appropriate sequence class for a given combination of
ALPHABET, AMBIGUITY and QUALITY."
  (let ((class (assocdr (list alphabet ambiguity quality)
                        *sequence-class-table* :test #'equal)))
    (unless class
      (error "invalid combination (~a ~a ~a)"
             alphabet ambiguity quality))
    class))

(defun make-seq (&rest initargs &key (alphabet :dna) ambiguity
                 &allow-other-keys)
  "Convenience constructor for bio-sequences."
  (let ((class (select-sequence-class alphabet ambiguity :default)))
    (multiple-value-bind (args remaining-initargs)
        (remove-args '(:alphabet :ambiguity) initargs)
      (apply #'make-instance class remaining-initargs))))

(defun make-quality-seq (&rest initargs &key (alphabet :dna) ambiguity
                         &allow-other-keys)
  "Convenience constructor for bio-sequences with quality."
  (let ((class (select-sequence-class alphabet ambiguity :quality)))
    (multiple-value-bind (args remaining-initargs)
        (remove-args '(:alphabet :ambiguity) initargs)
      (apply #'make-instance class remaining-initargs))))

(defmethod initialize-instance :after ((seq simple-dna-sequence) &key)
  (with-slots (token-seq length) seq
    (unless (encodedp token-seq '(unsigned-byte 2))
      (setf token-seq (encode-simple-seq token-seq #'encode-dna-2bit)))
    (setf length (length token-seq))))

(defmethod initialize-instance :after ((seq simple-rna-sequence) &key)
  (with-slots (token-seq length) seq
    (unless (encodedp token-seq '(unsigned-byte 2))
      (setf token-seq (encode-simple-seq token-seq #'encode-rna-2bit)))
    (setf length (length token-seq))))

(defmethod initialize-instance :after ((seq iupac-dna-sequence) &key)
  (with-slots (token-seq length) seq
    (unless (encodedp token-seq '(unsigned-byte 4))
      (setf token-seq (encode-iupac-seq token-seq #'encode-dna-4bit)))
    (setf length (length token-seq))))

(defmethod initialize-instance :after ((seq iupac-rna-sequence) &key)
  (with-slots (token-seq length) seq
    (unless (encodedp token-seq '(unsigned-byte 4))
      (setf token-seq (encode-iupac-seq token-seq #'encode-rna-4bit)))
    (setf length (length token-seq))))

(defmethod initialize-instance :after ((seq quality-mixin) &key)
  (with-slots (token-seq metric quality) seq
    (unless (= (length token-seq)
               (length quality))
      (error "token-seq and quality must be the same length but were (~a) and (~a) elements long"
             (length token-seq) (length quality)))
    (let ((decoder (cond ((eql :phred metric)
                          #'decode-phred-quality)
                         ((eql :illumina 
                               #'decode-illumina-quality))
                         (t
                          (error "invalid metric (~a), expected one of ~a"
                                 metric '(:phred :illumina))))))
      (setf quality (decode-quality quality decoder)))))

(defmethod length-of ((seq sequenced-mixin))
  (let ((token-seq (token-seq-of seq)))
    (length token-seq)))

(defmethod residue-of ((seq simple-dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-token seq index (encoded-tokens 2)))

(defmethod (setf residue-of) (value (seq simple-dna-sequence)
                              (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-token value seq index (encoded-tokens 2)))

(defmethod residue-of ((seq simple-rna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-token seq index (encoded-tokens 2)))

(defmethod (setf residue-of) (value (seq simple-rna-sequence)
                              (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-token value seq index (encoded-tokens 2)))

(defmethod residue-of ((seq iupac-dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-token seq index (encoded-tokens 4)))

(defmethod (setf residue-of) (value (seq iupac-dna-sequence)
                              (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-token value seq index (encoded-tokens 4)))

(defmethod residue-of ((seq iupac-rna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-token seq index (encoded-tokens 4)))

(defmethod (setf residue-of) (value (seq iupac-rna-sequence)
                              (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-token value seq index (encoded-tokens 4)))

(defmethod to-string ((seq simple-dna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (decode-token-array seq start end (encoded-tokens 2)))

(defmethod to-string ((seq simple-rna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (decode-token-array seq start end (encoded-tokens 2)))

(defmethod to-string ((seq iupac-dna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (decode-token-array seq start end (encoded-tokens 4)))

(defmethod to-string ((seq iupac-rna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (decode-token-array seq start end (encoded-tokens 4)))

(defmethod copy-sequence ((seq bio-sequence))
  (let ((token-seq (token-seq-of seq)))
    (make-instance (class-of seq) :token-seq
                   (make-array (length token-seq)
                               :element-type
                               (array-element-type token-seq)
                               :initial-contents token-seq))))

(defmethod print-object ((obj alphabet) stream)
  (format stream "<ALPHABET ~a>" (slot-value obj 'name)))

(defmethod print-object ((obj sequence-strand) stream)
  (with-slots (name token number) obj
      (format stream "<SEQUENCE-STRAND ~a/~a/~a>" name token number)))

(defmethod print-object ((obj simple-dna-sequence) stream)
  (print-seq-aux "SIMPLE-DNA-SEQUENCE" obj stream))

(defmethod print-object ((obj simple-rna-sequence) stream)
  (print-seq-aux "SIMPLE-RNA-SEQUENCE" obj stream))

(defmethod print-object ((obj iupac-dna-sequence) stream)
  (print-seq-aux "IUPAC-DNA-SEQUENCE" obj stream))

(defmethod print-object ((obj iupac-rna-sequence) stream)
  (print-seq-aux "IUPAC-RNA-SEQUENCE" obj stream))

(defmethod print-object ((obj simple-dna-quality-sequence) stream)
  (print-quality-seq-aux "SIMPLE-DNA-QUALITY-SEQUENCE" obj stream))

(defmethod print-object ((obj iupac-dna-quality-sequence) stream)
  (print-quality-seq-aux "IUPAC-DNA-QUALITY-SEQUENCE" obj stream))

(defun print-seq-aux (name obj stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of obj)))
    (if (and len (<= len *sequence-print-limit*))
        (format stream "<~a \"~a\">" name (to-string obj))
      (format stream "<~a length ~d>" name len))))

(defun print-quality-seq-aux (name obj stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of obj)))
    (if (and len (<= len *sequence-print-limit*))
        (format stream "<~a ~a quality, \"~a\">"
                name (metric-of obj) (to-string obj))
      (format stream "<~a ~a quality, length ~d>"
              name (metric-of obj) len))))

(defun encodedp (token-seq element-type)
  "Returns T if TOKEN-SEQ has elements encoded as ELEMENT-TYPE."
  (equal element-type (array-element-type token-seq)))
