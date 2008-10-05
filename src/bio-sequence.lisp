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

(deftype encoded-residues (n)
  `(simple-array (unsigned-byte ,n) *))

(deftype residues-length ()
  "Type for a sequence length."
  '(and fixnum (integer 1 *)))

(deftype quality-score ()
  "Type for sequence base quality score."
  '(signed-byte 8))

(defvar *sequence-print-limit* 50
  "Maximum length of sequence to be pretty-printed.")

(defmacro define-strand-decoder (type test forward reverse unknown)
  `(progn
    (defmethod decode-strand ((strand ,type) &key strict)
      (cond ((,test ,forward strand)
             *forward-strand*)
            ((,test ,reverse strand)
             *reverse-strand*)
            ((,test ,unknown strand)
             *unknown-strand*)
            (t
             (if strict
                 (error 'invalid-argument-error
                        :params 'strand
                        :args strand
                        :text (format nil "not a valid ~a strand designator"
                                      ',type))
               nil))))))

(define-strand-decoder string string= "+" "-" "?")
(define-strand-decoder fixnum = 1 -1 0)
(define-strand-decoder character char= #\+ #\- #\?)
(define-strand-decoder symbol eql :forward :reverse :unknown)

(defun find-alphabet (name)
  (multiple-value-bind (alphabet presentp)
      (gethash name *alphabets*)
    (unless presentp
      (error 'invalid-argument-error
             :params 'name
             :args name
             :text "no such alphabet"))
    alphabet))

(defun make-dna (residues &rest initargs
                 &key (encode t) &allow-other-keys)
  (let ((initargs (remove-args '(:encode) initargs)))
    (cond ((null residues)
           (apply #'make-instance 'virtual-dna-sequence initargs))
          (encode
           (make-encoded-vector-seq 'encoded-dna-sequence residues
                                    #'encode-dna-4bit '(unsigned-byte 4)
                                    initargs))
          (t
           (make-simple-vector-seq 'simple-dna-sequence residues
                                   initargs)))))

(defun make-rna (residues &rest initargs
                 &key (encode t) &allow-other-keys)
  (let ((initargs (remove-args '(:encode) initargs)))
    (cond ((null residues)
           (apply #'make-instance 'virtual-rna-sequence initargs))
          (encode
           (make-encoded-vector-seq 'encoded-rna-sequence residues
                                    #'encode-rna-4bit '(unsigned-byte 4)
                                    initargs))
          (t
           (make-simple-vector-seq 'simple-rna-sequence residues
                                   initargs)))))

(defun make-dna-quality (residues quality &rest initargs
                         &key (encode t) (metric :phred) &allow-other-keys)
  (unless (= (length residues)
             (length quality))
    (error 'invalid-argument-error
           :params '(residues quality)
           :args (list residues quality)
           :text "the residues and quality vectors were not the same length"))
  (let ((initargs (remove-args '(:encode :metric) initargs)))
    (cond (encode
           (apply #'make-instance 'dna-quality-sequence
                  :vector (ensure-encoded residues #'encode-dna-4bit
                                          '(unsigned-byte 4))
                  :quality (ensure-decoded-quality quality metric)
                  :metric metric
                  initargs))
          (t
            (error "Not implemented.")))))

(defun make-aa (residues &rest initargs
                &key (encode t) &allow-other-keys)
  (let ((initargs (remove-args '(:encode) initargs)))
    (cond ((null residues)
           (apply #'make-instance 'virtual-aa-sequence initargs))
          (encode
           (make-encoded-vector-seq 'encoded-aa-sequence residues
                                    #'encode-aa-7bit '(unsigned-byte 7)
                                    initargs))
          (t
           (make-simple-vector-seq 'simple-aa-sequence residues
                                   initargs)))))

(defun make-encoded-vector-seq (class residues encoder element-type initargs)
  (unless (and (vectorp residues) (not (zerop (length residues))))
    (error 'invalid-argument-error
           :params 'residues
           :args residues
           :text "expected a non-empty vector"))
  (apply #'make-instance class
         :vector (ensure-encoded residues encoder element-type) initargs))

(defun make-simple-vector-seq (class residues initargs)
  (declare (ignore class residues initargs))
  (error "Not implemented."))

;;; Printing methods
(defmethod print-object ((alphabet alphabet) stream)
  (format stream "#<ALPHABET ~a>" (slot-value alphabet 'name)))

(defmethod print-object ((strand sequence-strand) stream)
  (with-slots (name token number) strand
      (format stream "#<SEQUENCE-STRAND ~a/~a/~a>" name token number)))

(defmethod print-object ((seq dna-sequence) stream)
  (print-seq-aux "DNA-SEQUENCE" seq stream))

(defmethod print-object ((seq rna-sequence) stream)
  (print-seq-aux "RNA-SEQUENCE" seq stream))

(defmethod print-object ((seq dna-quality-sequence) stream)
  (print-quality-seq-aux "DNA-QUALITY-SEQUENCE" seq stream))

(defmethod print-object ((seq aa-sequence) stream)
  (print-seq-aux "AA-SEQUENCE" seq stream))

;;; Implementation methods
(defmethod anonymousp ((seq identity-mixin))
  (null (identity-of seq)))

(defmethod size-of ((alphabet alphabet))
  (length (tokens-of alphabet)))

(defmethod token-index ((alphabet alphabet) encoded-token)
  (gethash encoded-token (index-of alphabet)))

(defmethod memberp ((alphabet alphabet) (char character))
  (contains-char-p (tokens-of alphabet) char))

(defmethod explode-ambiguity ((alphabet (eql *dna*)) (char character))
  (explode-ambiguity-aux char #'encode-dna-4bit #'decode-dna-4bit))

(defmethod explode-ambiguity ((alphabet (eql *rna*)) (char character))
  (explode-ambiguity-aux char #'encode-rna-4bit #'decode-rna-4bit))

(defmethod ambiguousp ((alphabet alphabet) (char character))
  (> (length (explode-ambiguity alphabet char)) 1))

(defmethod simplep ((residues string) (alphabet alphabet))
  (loop
     for residue across residues
     never (ambiguousp alphabet residue)))

(defmethod strand-designator-p (strand)
  nil)

(defmethod strand-designator-p ((strand string))
  (decode-strand strand))

(defmethod strand-designator-p ((strand fixnum))
  (decode-strand strand))

(defmethod forward-strand-p ((strand (eql *forward-strand*)))
  t)

(defmethod reverse-strand-p ((strand (eql *reverse-strand*)))
  t)

(defmethod unknown-strand-p ((strand (eql *unknown-strand*)))
  t)

(defmethod invert-strand ((strand sequence-strand))
  (cond ((eql *forward-strand* strand)
         *reverse-strand*)
        ((eql *reverse-strand* strand)
         *forward-strand*)
        (t
         *unknown-strand*)))

(defmethod match-strand ((strand1 sequence-strand) (strand2 sequence-strand))
  (cond ((or (eql *unknown-strand* strand1) (eql *unknown-strand* strand2))
         t)
        (t
         (eql strand1 strand2))))

(defmethod virtualp ((seq token-sequence))
  (declare (ignore seq))
  nil)

(defmethod virtualp ((seq virtual-token-sequence))
  (declare (ignore seq))
  t)

(defmethod length-of ((seq vector-sequence))
  (with-slots (vector)
      seq
    (length vector)))

(defmethod single-stranded-p ((seq na-sequence))
  (with-slots (num-strands)
      seq
    (= 1 num-strands)))

(defmethod double-stranded-p ((seq na-sequence))
  (with-slots (num-strands)
      seq
    (= 2 num-strands)))

(defmethod residue-of ((seq encoded-dna-sequence) (index fixnum))
  (decode-dna-4bit (aref (vector-of seq) index)))

(defmethod (setf residue-of) (value (seq encoded-dna-sequence) (index fixnum))
  (with-slots (vector)
      seq
    (setf (aref vector index) (encode-dna-4bit value))))

(defmethod residue-of ((seq encoded-rna-sequence) (index fixnum))
  (with-slots (vector)
      seq
    (decode-rna-4bit (aref vector index))))

(defmethod (setf residue-of) (value (seq encoded-rna-sequence) (index fixnum))
  (with-slots (vector)
      seq
    (setf (aref vector index) (encode-rna-4bit value))))

(defmethod residue-of ((seq virtual-token-sequence) (index fixnum))
  (with-slots (length)
      seq
    (if (or (< index 0) (>= index length))
        (error 'invalid-argument-error
               :params 'index
               :args index
               :text "index must be >0 and < sequence length")
      *gap-char*)))

(defmethod num-gaps-of ((seq encoded-vector-sequence)
                        &key (start 0) (end (length-of seq)))
  (with-slots (vector)
      seq
    (loop
       for i from start below end
       count (= *encoded-gap-char* (aref vector i)))))

(defmethod to-string ((seq virtual-token-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (declare (ignore token-case))
  (make-string (- end start) :element-type 'base-char
               :initial-element *gap-char*))

(defmethod to-string ((seq encoded-dna-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (with-slots (vector)
      seq
    (let ((str (make-string (- end start) :element-type 'base-char))
          (seq-end (1- end))
          (str-start 0))
      (when (< 0 (length str))
        (copy-array vector start seq-end
                    str str-start #'decode-dna-4bit))
      (ecase token-case
        ((nil) str)
        (:lowercase str)
        (:uppercase (nstring-upcase str))))))

(defmethod to-string ((seq encoded-rna-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (with-slots (vector)
      seq
    (let ((str (make-string (- end start) :element-type 'base-char))
          (seq-end (1- end))
          (str-start 0))
      (when (< 0 (length str))
        (copy-array vector start seq-end
                    str str-start #'decode-rna-4bit))
      (ecase token-case
        ((nil) str)
        (:lowercase str)
        (:uppercase (nstring-upcase str))))))

(defmethod to-string ((seq encoded-aa-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (with-slots (vector)
      seq
    (let ((str (make-string (- end start) :element-type 'base-char))
          (seq-end (1- end))
          (str-start 0))
      (when (< 0 (length str))
        (copy-array vector start seq-end
                    str str-start #'decode-aa-7bit))
      (ecase token-case
        ((nil) str)
        (:lowercase str)
        (:uppercase (nstring-upcase str))))))

(defmethod subsequence ((seq vector-sequence) (start fixnum)
                        &optional end)
  (with-slots (vector)
      seq
    (make-instance (class-of seq)
                   :vector (token-subsequence vector start end))))

(defmethod subsequence  ((seq dna-quality-sequence) (start fixnum)
                         &optional end)
  (with-slots (vector quality metric)
      seq
    (make-instance 'dna-quality-sequence
                   :vector (token-subsequence vector start end)
                   :quality (subseq quality start end)
                   :metric metric)))

(defmethod subsequence ((seq virtual-token-sequence) (start fixnum)
                        &optional end)
  (with-slots (length)
      seq
    (let ((end (or end length)))
      (check-token-range length start end)
      (make-instance (class-of seq) :length (- end start)))))

(defmethod reverse-sequence ((seq vector-sequence))
  (with-slots (vector)
      seq 
    (make-instance (class-of seq) :vector (reverse vector))))

(defmethod reverse-sequence ((seq dna-quality-sequence))
  (with-slots (vector quality metric)
      seq
    (make-instance 'dna-quality-sequence
                   :vector (reverse vector)
                   :quality (reverse quality)
                   :metric metric)))

(defmethod reverse-sequence ((seq virtual-token-sequence))
  (with-slots (length)
      seq 
    (make-instance (class-of seq) :length length)))

(defmethod nreverse-sequence ((seq vector-sequence))
  (with-slots (vector)
      seq
    (setf vector (nreverse vector))
    seq))

(defmethod nreverse-sequence ((seq dna-quality-sequence))
  (with-slots (vector quality)
      seq
    (setf vector (nreverse vector)
          quality (nreverse quality))
    seq))

(defmethod nreverse-sequence ((seq virtual-token-sequence))
  seq)

(defmethod complement-sequence ((seq encoded-dna-sequence))
  (with-slots (vector)
      seq
    (make-instance 'encoded-dna-sequence :vector
                   (complement-tokens
                    (copy-seq vector) #'complement-dna-4bit))))

(defmethod complement-sequence ((seq encoded-rna-sequence))
  (with-slots (vector)
      seq
    (make-instance 'encoded-rna-sequence :vector
                   (complement-tokens
                    (copy-seq vector) #'complement-rna-4bit))))

(defmethod complement-sequence ((seq virtual-token-sequence))
  (with-slots (length)
      seq 
    (make-instance (class-of seq) :length length)))

(defmethod ncomplement-sequence ((seq encoded-dna-sequence))
  (with-slots (vector)
      seq
    (loop
       for i from 0 below (length vector)
       do (setf (aref vector i)
                (complement-dna-4bit (aref vector i)))))
  seq)

(defmethod ncomplement-sequence ((seq encoded-rna-sequence))
  (with-slots (vector)
      seq
    (loop
       for i from 0 below (length vector)
       do (setf (aref vector i)
                (complement-rna-4bit (aref vector i)))))
  seq)

(defmethod ncomplement-sequence ((seq virtual-token-sequence))
  seq)

(defmethod reverse-complement ((seq encoded-dna-sequence))
  (with-slots (vector)
      seq
    (make-instance 'encoded-dna-sequence :vector
                   (nreverse
                    (complement-tokens vector #'complement-dna-4bit)))))

(defmethod reverse-complement ((seq encoded-rna-sequence))
  (with-slots (vector)
      seq
    (make-instance 'encoded-rna-sequence :vector
                   (nreverse
                    (complement-tokens vector #'complement-rna-4bit)))))

(defmethod reverse-complement ((seq dna-quality-sequence))
  (with-slots (vector quality metric)
      seq
    (let ((s (make-instance 'dna-quality-sequence
                            :vector (copy-seq vector)
                            :quality (copy-seq quality)
                            :metric metric)))
      (nreverse-complement s))))

(defmethod reverse-complement ((seq virtual-token-sequence))
  (with-slots (length)
      seq 
    (make-instance (class-of seq) :length length)))

(defmethod nreverse-complement ((seq token-sequence))
  (nreverse-sequence (ncomplement-sequence seq)))

(defmethod nreverse-complement ((seq virtual-token-sequence))
  seq)

(defmethod search-sequence :around ((seq1 token-sequence)
                                    (seq2 token-sequence)
                                    &key from-end start1 start2 end1 end2)
  (declare (ignore from-end start1 start2 end1 end2))
  (if (subtypep (class-of (alphabet-of seq1))
                (class-of (alphabet-of seq2)))
      (call-next-method)
    nil))

(defmethod search-sequence ((seq1 vector-sequence)
                            (seq2 vector-sequence)
                            &key from-end start1 start2 end1 end2)
  (let ((vector1 (vector-of seq1))
        (vector2 (vector-of seq2))
        (start1 (or start1 0))
        (start2 (or start2 0)))
    (search vector1 vector2 :from-end from-end
                :start1 start1 :start2 start2
                :end1 end1 :end2 end2 :test #'eq)))

(defmethod residue-frequencies ((seq vector-sequence) (alphabet alphabet))
  (with-slots (vector)
      seq
    (let ((frequencies (make-array (length (tokens-of alphabet))
                                   :element-type 'fixnum :initial-element 0)))
      (loop
         for elt across vector
         do (incf (aref frequencies (token-index alphabet elt))))
      (pairlis (coerce (copy-seq (tokens-of alphabet)) 'list)
               (coerce frequencies 'list)))))

;;; Utility functions
(defun print-seq-aux (name seq stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of seq)))
    (if (<= len *sequence-print-limit*)
        (format stream "#<~a \"~a\">" name (to-string seq))
      (format stream "#<~a length ~d>" name len))))

(defun print-quality-seq-aux (name seq stream)
  "Helper function for printing bio-sequence objects."
  (with-slots (quality metric)
      seq
    (let ((len (length-of seq)))
      (if (<= len *sequence-print-limit*)
          (format stream "#<~a \"~a\" ~a quality \"~a\">"
                  name (to-string seq) metric
                  (quality-string quality metric))
        (format stream "#<~a ~a quality, length ~d>"
                name metric len)))))

(defun quality-string (quality metric)
  "Wrapper for ENCODE-QUALITY that encodes QUALITY, an array of bytes
representing base quality scores, as a string using an encoder
appropriate for METRIC, a quality metric \(:PHRED or :ILLUMINA\)."
  (let ((encoder (ecase metric
                   (:phred #'encode-phred-quality)
                   (:illumina #'encode-illumina-quality))))
    (encode-quality quality encoder)))

(defun encode-quality (quality encoder)
  "Encodes QUALITY, an array of bytes representing base quality
scores, as a string using function ENCODER."
  (let ((quality-str (make-string (length quality)
                                  :element-type 'base-char)))
    (copy-array quality 0 (1- (length quality))
                quality-str 0 encoder)
    quality-str))

(defun decode-quality (quality decoder)
  "Decodes the QUALITY, a string, as into a new array using function
DECODER."
  (let ((quality-seq (make-array (length quality)
                                 :element-type 'quality-score)))
    (copy-array quality 0 (1- (length quality))
                quality-seq 0 decoder)
    quality-seq))

(defun ensure-encoded (vector encoder element-type)
  "If VECTOR is a string, returns an encoded token vector of element
type \(unsigned-byte 4\) of the same length, otherwise returns
VECTOR. ENCODER is the encoding function used to convert characters to
\(unsigned-byte 4\)."
  (if (stringp vector)
      (let ((encoded (make-array (length vector)
                                 :element-type element-type)))
        (when (< 0 (length vector))
          (copy-array vector 0 (1- (length encoded))
                      encoded 0 encoder))
        encoded)
    vector))

(defun ensure-decoded-quality (quality metric)
  "If QUALITY is a string, returns a decoded quality vector of
integers of the same length, otherwise returns QUALITY. METRIC is
either :PHRED or :ILLUMINA, denoting the quality metric to be used
when decoding."
  (if (stringp quality)
      (let ((decoder (ecase metric
                       (:phred #'decode-phred-quality)
                       (:illumina #'decode-illumina-quality))))
        (decode-quality quality decoder))
    quality))

(defun token-subsequence (tokens start end)
  "Returns a subsequence of TOKENS between indices START and END."
  (let ((end (or end (length tokens))))
    (check-token-range (length tokens) start end)
    (let ((sub-seq (make-array (- end start)
                               :element-type
                               (array-element-type tokens))))
      (copy-array tokens start (1- end)
                  sub-seq 0)
      sub-seq)))

(defun check-token-range (length start end)
  (cond ((or (< start 0) (> start length))
         (error 'invalid-argument-error
                :params 'start
                :args start
                :text "start index must be >0 and <= sequence length"))
        ((or (< end 0) (> end length))
         (error 'invalid-argument-error
                :params 'end
                :args end
                :text "end index must be >0 and <= sequence length"))
        ((< end start)
         (error 'invalid-argument-error
                :params '(start end)
                :args (list start end)
                :text "end index must be equal to or greater than start index"))
        (t t)))

(defun complement-tokens (tokens comp-fn &optional (start 0) end)
  "Returns a complemented copy of RESIDUES populated with elements
from RESIDUES that have been transformed by COMP-FN, starting at the
first element, or index START, and continuing to the last residue, or
index END."
  (let* ((end (or end (length tokens)))
         (comp-seq (make-array (- end start)
                               :element-type
                               (array-element-type tokens))))
    (copy-array tokens start (1- end)
                comp-seq 0 comp-fn)
    comp-seq))

(defun explode-ambiguity-aux (char encoder decoder)
  "Returns a list of the ambiguity characters represented by CHAR."
  (let ((encoded (funcall encoder char)))
   (loop
      for b from 0 below (integer-length encoded)
      when (logbitp b encoded)
      collect (funcall decoder (ash 1 b)) into exploded
      finally (return (sort exploded #'char<=)))))
