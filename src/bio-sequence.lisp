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

(deftype encoded-tokens (n)
  `(simple-array (unsigned-byte ,n) *))

(deftype token-seq-length ()
  "Type for a sequence length."
  '(and fixnum (integer 1 *)))

(defvar *sequence-class-table*
  (pairlis '((:dna nil :default) (:rna nil :default)
             (:dna :iupac :default) (:rna :iupac :default)
             (:dna nil :quality)
             (:dna :iupac :quality))
           '(simple-dna-sequence simple-rna-sequence
             iupac-dna-sequence iupac-rna-sequence
             simple-dna-quality-sequence
             iupac-dna-quality-sequence))
  "Mappings of (alphabet ambuiguity sequence-quality) tuples to CLOS
bio-sequence classes.")

(defvar *sequence-print-limit* 50
  "Maximum length of sequence to be pretty-printed.")


(defmacro encode-token (value seq index token-seq-type)
  "Sets the residue token VALUE at INDEX in SEQ. The residue tokens
are encoded as an array of type TOKEN-SEQ-TYPE in SEQ."
  (let ((token-seq (gensym))
        (encoder (gensym)))
    `(let ((,token-seq (token-seq-of ,seq))
           (,encoder (encoder-of (alphabet-of ,seq))))
       (declare (type ,token-seq-type ,token-seq)
                (type function ,encoder))
       (setf (aref ,token-seq ,index)
             (funcall ,encoder ,value)))))

(defun encode-token-2bit (value seq index)
  "Sets the residue token VALUE at INDEX in SEQ. The residue tokens
are encoded as an array of type (encoded-tokens 2) in SEQ."
  (let ((token-seq (token-seq-of seq))
        (encoder (encoder-of (alphabet-of seq))))
    (declare (type (encoded-tokens 2) token-seq)
             (type function encoder))
    (setf (aref token-seq index)
          (funcall encoder value))))

(defmacro decode-token (seq index token-seq-type)
  "Returns the decoded residue token from INDEX in SEQ. The residue
tokens are encoded as an array of type TOKEN-SEQ-TYPE in SEQ."
  (let ((token-seq (gensym))
        (decoder (gensym)))
    `(let ((,token-seq (token-seq-of ,seq))
           (,decoder (decoder-of (alphabet-of seq))))
       (declare (type ,token-seq-type ,token-seq)
                (type function ,decoder))
       (funcall ,decoder (aref ,token-seq ,index)))))

(defmacro decode-token-array (seq start end token-seq-type)
  "Returns a simple-base-string representing the token-seq of SEQ
from residues START to END, inclusive. The residues are encoded as
type TOKEN-SEQ-TYPE in the token-seq slot of SEQ."
  (let ((token-seq (gensym))
        (decoder (gensym))
        (dest (gensym))
        (source-end (gensym))
        (dest-start (gensym)))
    `(let ((,token-seq (token-seq-of ,seq))
           (,decoder (decoder-of (alphabet-of seq)))
           (,dest (make-string (- ,end ,start)
                               :element-type 'base-char))
           (,source-end (1- ,end))
           (,dest-start 0))
      (declare (type ,token-seq-type ,token-seq)
               (type function ,decoder)
               (type simple-base-string ,dest)
               (type fixnum ,source-end ,dest-start))
      (when (< 0 (length ,dest))
        (copy-array ,token-seq ,start ,source-end
                    ,dest ,dest-start ,decoder))
      ,dest)))

(defun encode-2bit (str encoder)
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

(defun encode-4bit (str encoder)
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

(defun ensure-2bit-encoded (vector encoder)
  "If VECTOR is not of element-type (unsigned-byte 2), attempts to
encode it as such with ENCODER."
  (if (equal (array-element-type vector) '(unsigned-byte 2))
      vector
    (encode-2bit vector encoder)))

(defun ensure-4bit-encoded (vector encoder)
  "If VECTOR is not of element-type (unsigned-byte 4), attempts to
encode it as such with ENCODER."
  (if (equal (array-element-type vector) '(unsigned-byte 4))
      vector
    (encode-4bit vector encoder)))

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
      (error (msg "Invalid (alphabet ambiguity quality)"
                  "combination (~a ~a ~a).")
             alphabet ambiguity quality))
    class))

(defun make-seq (&rest initargs &key (alphabet :dna) ambiguity
                 &allow-other-keys)
  "
 Purpose: Convenience constructor for bio-sequences. 

     key: ALPHABET (:dna or :rna).
          AMBIGUITY (:iupac or NIL).
   &rest: INITARGS passed to MAKE-INSTANCE.

 Returns: A BIO-SEQUENCE object.
"
  (let ((class (select-sequence-class alphabet ambiguity :default)))
    (multiple-value-bind (args remaining-initargs)
        (remove-args '(:alphabet :ambiguity) initargs)
      (declare (ignore args))
      (apply #'make-instance class remaining-initargs))))

(defun make-quality-seq (&rest initargs &key (alphabet :dna) ambiguity
                         &allow-other-keys)
  "
 Purpose: Convenience constructor for bio-sequences with quality. 

     key: ALPHABET (:dna or :rna).
          AMBIGUITY (:iupac or NIL).
   &rest: INITARGS passed to MAKE-INSTANCE.

 Returns: A BIO-SEQUENCE object.
"
  (let ((class (select-sequence-class alphabet ambiguity :quality)))
    (multiple-value-bind (args remaining-initargs)
        (remove-args '(:alphabet :ambiguity) initargs)
      (declare (ignore args))
      (apply #'make-instance class remaining-initargs))))

(defmethod anonymousp ((obj identity-mixin))
  (null (identity-of obj)))

(defmethod simplep ((token-seq string) (alphabet (eql :dna)))
  (let ((simple (tokens-of (find-alphabet :dna))))
    (loop
       for token across token-seq
       always (find (char-downcase token) simple))))

(defmethod simplep ((token-seq string) (alphabet (eql :rna)))
  (let ((simple (tokens-of (find-alphabet :rna))))
    (loop
       for token across token-seq
       always (find (char-downcase token) simple))))

(defmethod slot-unbound (class (obj alphabet) (slot (eql 'encoded-index)))
  (let ((index-table (make-hash-table))
        (encoder (encoder-of obj))
        (tokens (tokens-of obj)))
    (loop
       for i from 0 below (length tokens)
       do (setf (gethash (funcall encoder (aref tokens i)) index-table) i))
    (setf (slot-value obj 'encoded-index)
          (lambda (element)
            (gethash element index-table)))))

(defmethod initialize-instance :after ((seq simple-dna-sequence) &key)
  (initialize-seq seq #'ensure-2bit-encoded (encoder-of (alphabet-of seq))))

(defmethod initialize-instance :after ((seq simple-rna-sequence) &key)
  (initialize-seq seq #'ensure-2bit-encoded (encoder-of (alphabet-of seq))))

(defmethod initialize-instance :after ((seq iupac-dna-sequence) &key)
  (initialize-seq seq #'ensure-4bit-encoded (encoder-of (alphabet-of seq))))

(defmethod initialize-instance :after ((seq iupac-rna-sequence) &key)
  (initialize-seq seq #'ensure-4bit-encoded (encoder-of (alphabet-of seq))))

(defmethod initialize-instance :after ((seq quality-mixin) &key)
  (with-slots (token-seq metric quality) seq
    (unless (= (length token-seq)
               (length quality))
      (error (msg "Token-seq and quality must be the same length"
                  "but were ~a and ~a elements long, respectively.")
             (length token-seq) (length quality)))
    (let ((decoder (cond ((eql :phred metric)
                          #'decode-phred-quality)
                         ((eql :illumina
                               #'decode-illumina-quality))
                         (t
                          (error "Invalid metric ~a: expected one of ~a."
                                 metric '(:phred :illumina))))))
      (setf quality (decode-quality quality decoder)))))

(defmethod token-seq-of ((seq bio-sequence))
  (slot-value seq 'token-seq))

(defmethod length-of ((seq bio-sequence))
  (with-slots (token-seq length) seq
    (if (null token-seq)
        length
      (length token-seq))))

(defmethod (setf length-of) (value (seq bio-sequence))
  (with-slots (token-seq length) seq
    (if (null token-seq)
        (setf length value)
      (error (msg "Invalid operation: the length of a concrete sequence"
                  "may not be changed.")))))

(defmethod virtualp ((seq bio-sequence))
  (null (slot-value seq 'token-seq)))

(defmethod residue-of :around ((seq bio-sequence) index)
  (when (virtualp seq)
    (error (msg "Invalid operation: cannot access a residue"
                "in a virtual sequence.")))
  (call-next-method))

(defmethod (setf residue-of) :around (value (seq bio-sequence) index)
  (when (virtualp seq)
    (error (msg "Invalid operation: cannot access a residue"
                "in a virtual sequence.")))
  (call-next-method))

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

;; (defun cumulative-lengths (ranges)
;;   (loop
;;      for range in ranges
;;      for n = (length-of range) then (+ n (length-of range))
;;      collect n))

;; (defun map-position (source-position ranges-in-source ranges-in-target)
;;   "Assumes all ranges are congruent."
;;   (let* ((target-ranges (sort (copy-seq ranges-in-target)
;;                               #'< :key #'start-of))
;;          (cumul-lengths (cumulative-lengths target-ranges))
;;          (source-ranges (loop
;;                            for len in cumul-lengths
;;                            and start = 0 then len
;;                            collect (cons start len)))
;;          (range-num (position source-position cumul-lengths :test #'<)))
;;     (+ (start-of (elt target-ranges range-num))
;;        (- source-position (start-of (elt source-ranges range-num))))))


(defmethod to-string :around ((seq bio-sequence) &optional start end)
  (declare (ignore start end))
  (if (virtualp seq)
      "<virtual>"
    (call-next-method)))

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

;;; FIXME -- add a means of making a reversed and/or complemented view
;;; of a sequence without modifying it

(defun reverse-complement-index (index len)
  (- len index))

(defmethod subsequence ((seq bio-sequence) (start fixnum) &optional end)
  (let* ((token-seq (token-seq-of seq))
         (end (or end (length token-seq)))
         (sub-seq (make-array (- end start)
                              :element-type
                              (array-element-type token-seq))))
    (copy-array token-seq start (1- end)
                sub-seq 0)
    (make-instance (class-of seq) :token-seq sub-seq)))

(defmethod reverse-sequence ((seq bio-sequence))
  (make-instance (class-of seq)
                 :token-seq (reverse (token-seq-of seq))))

(defmethod nreverse-sequence ((seq bio-sequence))
  (make-instance (class-of seq)
                 :token-seq (nreverse (token-seq-of seq))))

;;; FIXME -- factor out the common code in the complement methods,
;;; perhaps when or if it's time to add type declarations
(defun complement-token-seq (token-seq comp-fn &optional (start 0) end)
  (let* ((end (or end (length token-seq)))
         (comp-seq (make-array (- end start)
                               :element-type
                               (array-element-type token-seq))))
    (copy-array token-seq start (1- end)
                comp-seq 0 comp-fn)
    comp-seq))

(defmethod complement-sequence ((seq simple-dna-sequence)
                                &optional (start 0) end)
  (make-instance 'simple-dna-sequence :token-seq
                 (complement-token-seq
                  (token-seq-of seq) #'complement-dna-2bit start end)))

(defmethod complement-sequence ((seq simple-dna-sequence)
                                &optional (start 0) end)
  (make-instance 'iupac-dna-sequence :token-seq
                 (complement-token-seq
                  (token-seq-of seq) #'complement-dna-4bit start end)))

(defmethod reverse-complement ((seq simple-dna-sequence)
                               &optional (start 0) end)
  (make-instance 'simple-dna-sequence :token-seq
                 (nreverse
                  (complement-token-seq
                   (token-seq-of seq) #'complement-dna-2bit start end))))

(defmethod reverse-complement ((seq iupac-dna-sequence)
                               &optional (start 0) end)
  (make-instance 'iupac-dna-sequence :token-seq
                 (nreverse
                  (complement-token-seq
                   (token-seq-of seq) #'complement-dna-4bit start end))))

(defmethod residue-frequencies :before ((seq bio-sequence))
  (when (virtualp seq)
    (error "Invalid operation on virtual sequence ~a." seq)))

(defmethod residue-frequencies ((seq bio-sequence))
  (let ((index-fn (encoded-index-of (alphabet-of seq)))
        (frequencies (make-array (length (tokens-of (alphabet-of seq)))
                                 :element-type 'fixnum :initial-element 0))
        (token-seq (token-seq-of seq)))
    (loop
       for token across token-seq
       do (incf (aref frequencies (funcall index-fn token))))
    (pairlis (coerce (copy-seq (tokens-of (alphabet-of seq))) 'list)
             (coerce frequencies 'list))))


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
    (if (and len (<= len *sequence-print-limit*)
             (not (virtualp obj)))
        (format stream "<~a \"~a\">" name (to-string obj))
      (format stream "<~a length ~d>" name len))))

(defun print-quality-seq-aux (name obj stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of obj)))
    (if (and len (<= len *sequence-print-limit*)
             (not (virtualp obj)))
        (format stream "<~a ~a quality, \"~a\">"
                name (metric-of obj) (to-string obj))
      (format stream "<~a ~a quality, length ~d>"
              name (metric-of obj) len))))

(defun process-token-seq-args (token-seq length)
  "Returns its arguments, having checked their consistency for use
when making bio-sequence instances."
  (cond ((and (null token-seq)
              (null length))
         (error "Invalid token-seq and length: expected one to be non-NIL."))
        ((and token-seq
              (null length))
         (unless (and (vectorp token-seq)
                      (not (zerop (length token-seq))))
           (error "Invalid token-seq: expected a non-empty vector."))
         (values token-seq (length token-seq)))
        ((and (null token-seq)
              length)
         (unless (typep length 'token-seq-length)
           (error "Invalid length: expected a fixnum >= 1."))
         (values token-seq length))
        (t
         (error "Invalid token-seq and length: expected one to be NIL."))))

(defun initialize-seq (seq seq-encoder token-encoder)
  "Returns SEQ, having initialized the token-seq of bio-sequence SEQ
using SEQ-ENCODER to encode the vector and TOKEN-ENCODER to encode the
elements therein."
  (with-slots (token-seq length) seq
    (multiple-value-bind (valid-token-seq valid-length)
        (process-token-seq-args token-seq length)
      (if valid-token-seq
          (setf token-seq (funcall seq-encoder valid-token-seq token-encoder))
        (setf length valid-length))))
  seq)
