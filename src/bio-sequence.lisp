
;;; Glossary
;;;
;;; Residue token: a character symbolising a biological sequence residue
;;; in standard nomenclature e.g. t, c, a and g for DNA.
;;;
;;; Encoded token: a residue token encoded in some numeric format e.g.
;;; a 2-bit byte for DNA.
;;;

(in-package :bio-sequence)

(defvar *seq-len-print-limit* 50)

(deftype encoded-tokens (n)
  `(simple-array (unsigned-byte ,n) *))

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
       (gpu:copy-array ,token-seq ,start ,source-end
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
    (gpu:copy-array str 0 (1- (length str))
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
    (gpu:copy-array str 0 (1- (length str))
                    token-seq 0 encoder)
    token-seq))

(defun decode-quality (quality decoder)
  (let ((quality-seq (make-array (length quality)
                                 :element-type '(unsigned-byte 8))))
    (gpu:copy-array quality 0 (1- (length quality))
                    quality-seq 0 decoder)
    quality-seq))

(defun make-simple-seq (class str encoder &rest initargs)
  "Returns a new bio-sequence of CLASS composed of the tokens in the
simple-string STR encoded as (unsigned-byte 2) with ENCODER."
  (apply #'make-instance class
         :token-seq (encode-simple-seq str encoder) initargs))

(defun make-iupac-seq (class str encoder &rest initargs)
  "Returns a new bio-sequence of CLASS composed of the tokens in the
simple-string STR encoded as (unsigned-byte 4) with ENCODER."
  (apply #'make-instance class
         :token-seq (encode-iupac-seq str encoder) initargs))

(defun make-dna-seq (str &key ambiguity)
   "Returns a new DNA-SEQUENCE object with residues specified by
simple-string STR. Base ambiguity may be defined with the :AMBIGUITY
key. Accepted values for :AMBGUITY are NIL (no ambiguity, the default)
and :IUPAC (IUPAC ambiguity)."
   (cond ((null ambiguity)
          (make-simple-seq 'simple-dna-sequence
                           str #'encode-dna-2bit))
         ((eql :iupac ambiguity)
          (make-iupac-seq 'iupac-dna-sequence
                          str #'encode-dna-4bit))
         (t
          (error "Illegal ambiguity: ~a. Expected one of ~a"
                 ambiguity '(nil :iupac)))))

(defun make-rna-seq (str &key ambiguity)
   "Returns a new RNA-SEQUENCE object with residues specified by
simple-string STR. Base ambiguity may be defined with the :AMBIGUITY
key. Accepted values for :AMBGUITY are NIL (no ambiguity, the default)
and :IUPAC (IUPAC ambiguity)."
   (cond ((null ambiguity)
          (make-simple-seq 'simple-rna-sequence
                           str #'encode-rna-2bit))
         ((eql :iupac ambiguity)
          (make-iupac-seq 'iupac-rna-sequence
                          str #'encode-rna-4bit))
         (t
          (error "Illegal ambiguity: ~a. Expected one of ~a"
                 ambiguity '(nil :iupac)))))

(defun make-dna-quality-seq (str quality
                             &key ambiguity (metric :phred))
  "Returns a new DNA-SEQUENCE object with residues specified by
simple-string STR and base quality specified by simple-string
QUALITY. Base ambiguity may be defined with the :AMBIGUITY
key. Accepted values for :AMBGUITY are NIL (no ambiguity, the default)
and :IUPAC (IUPAC ambiguity). The quality metric may be defined with
the :METRIC key. Accepted values for :METRIC are :PHRED (Phred
quality, the default) or :ILLUMINA (Illumina quality)."
  (let ((qual-decoder
         (cond ((eql :phred metric)
                #'decode-phred-quality)
               ((eql :illumina metric)
                #'decode-illumina-quality)
               (t
                (error "Illegal metric: ~a. Expected one of ~a"
                       metric '(:phred :illumina))))))
    (cond ((null ambiguity)
           (make-simple-seq 'simple-dna-quality-sequence
                            str #'encode-dna-2bit
                            :metric metric
                            :quality (decode-quality quality
                                                     qual-decoder)))
          ((eql :iupac ambiguity)
           (make-iupac-seq 'iupac-dna-quality-sequence
                           str #'encode-dna-4bit
                           :metric metric
                           :quality (decode-quality quality
                                                    qual-decoder)))
          (t
           (error "Illegal ambiguity: ~a. Expected one of ~a"
                  metric '(nil :iupac))))))

(defmethod length-of ((seq bio-sequence))
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
    (if (<= len *seq-len-print-limit*)
        (format stream "<~a \"~a\">" name (to-string obj))
      (format stream "<~a length ~d>" name len))))

(defun print-quality-seq-aux (name obj stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of obj)))
    (if (<= len *seq-len-print-limit*)
        (format stream "<~a ~a quality, \"~a\">"
                name (metric-of obj) (to-string obj))
      (format stream "<~a ~a quality, length ~d>"
              name (metric-of obj) len))))
