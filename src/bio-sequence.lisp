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

(deftype quality-score ()
  "Type for sequence base quality score."
  '(signed-byte 8))

(defparameter *sequence-print-limit* 160
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

(defun make-dna (residues &rest initargs
                 &key (encode t) &allow-other-keys)
  "Returns a new DNA sequence object.

Arguments:

- residues (vector or NIL): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.

Rest:

- initargs: Any initialization arguments.

Key:

- encode (boolean): Indicates whether the residues should be
  encoded.

Returns:

- A DNA sequence object."
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
  "Returns a new RNA sequence object.

Arguments:

- residues (vector or NIL): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.

Rest:

- initargs: Any initialization arguments.

Key:

- encode (boolean): Indicates whether the residues should be
  encoded.

Returns:

- A RNA sequence object."
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
  "Returns a new DNA sequence object with quality.

Arguments:

- residues (vector or NIL): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.
- quality (vector integer): A vector of integers, the same length as
  RESIDUES, containing quality data.

Rest:

- initargs: Any initialization arguments.

Key:

- encode (boolean): Indicates whether the residues should be
  encoded.
- metric (symbol): A symbol indicating the quality metric used in
  QUALITY. Either :phred or :illumina.

Returns:

- A DNA sequence object."
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
  "Returns a new AA sequence object.

Arguments:

- residues (vector or NIL): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.

Rest:

- initargs: Any initialization arguments.

Key:

- encode (boolean): Indicates whether the residues should be
  encoded.

Returns:

- An AA sequence object."
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
  (unless (vectorp residues)
    (error 'invalid-argument-error
           :params 'residues
           :args residues
           :text "expected a vector"))
  (apply #'make-instance class
         :vector (ensure-encoded residues encoder element-type) initargs))

(defun make-simple-vector-seq (class residues initargs)
  (declare (ignore class residues initargs))
  (error "Not implemented."))

(defun bio-sequence-p (obj)
  "Returns T if OBJ is a subtype of {defclass bio-sequence} , or NIL
otherwise."
  (subtypep (class-of obj) 'bio-sequence))

(defun na-sequence-p (obj)
   "Returns T if OBJ is a subtype of {defclass na-sequence} , or NIL
otherwise."
  (subtypep (class-of obj) 'na-sequence))

(defun dna-sequence-p (obj)
  "Returns T if OBJ is a subtype of {defclass dna-sequence} , or NIL
otherwise."
  (subtypep (class-of obj) 'dna-sequence))

(defun rna-sequence-p (obj)
  "Returns T if OBJ is a subtype of {defclass rna-sequence} , or NIL
otherwise."
  (subtypep (class-of obj) 'rna-sequence))

(defun aa-sequence-p (obj)
  "Returns T if OBJ is a subtype of {defclass aa-sequence} , or NIL
otherwise."
  (subtypep (class-of obj) 'aa-sequence))

(defun same-biotype-p (&rest seqs)
  "Returns T if SEQS are {defclass bio-sequence} s that share exactly
the same alphabet, or NIL otherwise."
  (cond ((null seqs)
         nil)
        ((not (bio-sequence-p (first seqs)))
         nil)
        (t
         (loop
            with alphabet = (alphabet-of (first seqs))
            for seq in (rest seqs)
            always (and (bio-sequence-p seq)
                        (eql alphabet (alphabet-of seq)))))))

(defun same-strand-num-p (&rest seqs)
  "Returns T if SEQS are {defclass na-sequence} s that share the same
number of strands, or NIL otherwise."
  (if (null seqs)
         nil
    (loop
       with num-strands = (num-strands-of (first seqs))
       for seq in (rest seqs)
       always (= num-strands (num-strands-of seq)))))

(defun concat-sequence (&rest seqs)
  (if (apply #'same-biotype-p seqs)
      (let ((construct (cond ((dna-sequence-p (first seqs))
                              #'make-dna)
                             ((rna-sequence-p (first seqs))
                              #'make-rna)
                             ((aa-sequence-p (first seqs))
                              #'make-aa)
                             (t
                              (error 'invalid-argument-error
                                     :params 'seqs
                                     :args seqs
                                     :text (msg "expected all sequences to be"
                                                "one of DNA, RNA or AA"))))))
        ;; We use coerce-sequence to avoid a special case for virtual
        ;; sequences. We can't simply concatenate the encoded residue
        ;; vectors because virtual sequences do not have them
        (funcall construct (apply #'concatenate 'string
                                  (mapcar (lambda (s)
                                            (coerce-sequence s 'string))
                                          seqs))))
    (error 'invalid-argument-error
           :params 'seqs
           :args seqs
           :text "expected all sequences to be one of DNA, RNA or AA")))

;;; Printing methods
(defmethod print-object ((alphabet alphabet) stream)
  (format stream "#<ALPHABET ~a>" (slot-value alphabet 'name)))

(defmethod print-object ((strand sequence-strand) stream)
  (with-accessors ((name name-of) (token token-of) (number number-of))
      strand
    (format stream "#<SEQUENCE-STRAND ~a/~a/~a>" name token number)))

(defmethod print-object ((seq dna-sequence) stream)
  (%print-seq "DNA-SEQUENCE" seq stream))

(defmethod print-object ((seq rna-sequence) stream)
  (%print-seq "RNA-SEQUENCE" seq stream))

(defmethod print-object ((seq dna-quality-sequence) stream)
  (%print-quality-seq "DNA-QUALITY-SEQUENCE" seq stream))

(defmethod print-object ((seq aa-sequence) stream)
  (%print-seq "AA-SEQUENCE" seq stream))

;;; Implementation methods
(defmethod anonymousp ((seq identity-mixin))
  (null (identity-of seq)))

(defmethod size-of ((alphabet alphabet))
  (length (tokens-of alphabet)))

(defmethod token-index ((alphabet alphabet) encoded-token)
  (gethash encoded-token (index-of alphabet)))

(defmethod memberp ((alphabet alphabet) (char character))
  (find char (tokens-of alphabet) :test #'char=))

(defmethod explode-ambiguity ((alphabet (eql *dna*)) (char character))
  (%explode-ambiguity char #'encode-dna-4bit #'decode-dna-4bit))

(defmethod explode-ambiguity ((alphabet (eql *rna*)) (char character))
  (%explode-ambiguity char #'encode-rna-4bit #'decode-rna-4bit))

(defmethod strand-designator-p (strand)
  nil)

(defmethod strand-designator-p ((strand string))
  (decode-strand strand))

(defmethod strand-designator-p ((strand fixnum))
  (decode-strand strand))

(defmethod strand-designator-p ((strand character))
  (decode-strand strand))

(defmethod strand-designator-p ((strand symbol))
  (decode-strand strand))

(defmethod forward-strand-p ((strand (eql *forward-strand*)))
  t)

(defmethod forward-strand-p ((strand (eql *reverse-strand*)))
  nil)

(defmethod forward-strand-p ((strand (eql *unknown-strand*)))
  nil)

(defmethod reverse-strand-p ((strand (eql *forward-strand*)))
  nil)

(defmethod reverse-strand-p ((strand (eql *reverse-strand*)))
  t)

(defmethod reverse-strand-p ((strand (eql *unknown-strand*)))
  nil)

(defmethod unknown-strand-p ((strand (eql *forward-strand*)))
  nil)

(defmethod unknown-strand-p ((strand (eql *reverse-strand*)))
  nil)

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

(defmethod ambiguousp ((seq virtual-token-sequence))
  t)

(defmethod ambiguousp ((seq vector-sequence))
  (with-accessors ((alphabet alphabet-of))
      seq
    (%ambiguousp seq alphabet)))

(defmethod simplep ((seq virtual-token-sequence))
  nil)

(defmethod simplep ((seq vector-sequence))
  (with-accessors ((alphabet alphabet-of))
      seq
    (not (%ambiguousp seq alphabet))))

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

(defmethod num-strands-of ((seq aa-sequence))
  (error 'bio-sequence-op-error
         :text "Amino-acid sequences are not stranded."))

(defmethod element-of ((seq encoded-dna-sequence) (index fixnum))
  (decode-dna-4bit (aref (vector-of seq) index)))

(defmethod (setf element-of) (value (seq encoded-dna-sequence) (index fixnum))
  (with-slots (vector)
      seq
    (setf (aref vector index) (encode-dna-4bit value))))

(defmethod element-of ((seq encoded-rna-sequence) (index fixnum))
  (decode-rna-4bit (aref (vector-of seq) index)))

(defmethod (setf element-of) (value (seq encoded-rna-sequence) (index fixnum))
  (with-accessors ((vector vector-of))
      seq
    (setf (aref vector index) (encode-rna-4bit value))))

(defmethod element-of ((seq virtual-dna-sequence) (index fixnum))
  (with-accessors ((length length-of))
      seq
    (%check-token-range length index index)
    #\n))

(defmethod element-of ((seq virtual-rna-sequence) (index fixnum))
  (with-accessors ((length length-of))
      seq
    (%check-token-range length index index)
    #\n))

(defmethod element-of ((seq virtual-aa-sequence) (index fixnum))
  (with-accessors ((length length-of))
      seq
    (%check-token-range length index index)
    #\X))

(defun residue-of (seq index)
  "Returns the residue of SEQ at INDEX. This is a synonym of ELEMENT-OF."
  (element-of seq index))

(defun (setf residue-of) (value seq index)
  "Sets the residue of SEQ at INDEX to VALUE. This is a synonym of
\(SETF ELEMENT-OF\)."
  (setf (element-of seq index) value))

(defmethod num-gaps-of ((seq encoded-vector-sequence)
                        &key (start 0) (end (length-of seq)))
  (with-accessors ((vector vector-of))
      seq
    (loop
       for i from start below end
       count (= +encoded-gap-char+ (aref vector i)))))

(defmethod coerce-sequence ((seq virtual-dna-sequence) (type (eql 'string))
                            &key (start 0) (end (length-of seq)))
  (%to-string-virtual seq #\n start end :lowercase))

(defmethod coerce-sequence ((seq virtual-rna-sequence) (type (eql 'string))
                            &key (start 0) (end (length-of seq)))
  (%to-string-virtual seq #\n start end :lowercase))

(defmethod coerce-sequence ((seq virtual-aa-sequence) (type (eql 'string))
                            &key (start 0) (end (length-of seq)))
  (%to-string-virtual seq #\X start end :uppercase))

(defmethod coerce-sequence ((seq encoded-dna-sequence) (type (eql 'string))
                            &key (start 0) (end (length-of seq)))
  (%to-string-encoded seq #'decode-dna-4bit start end :lowercase))

(defmethod coerce-sequence ((seq encoded-rna-sequence) (type (eql 'string))
                            &key (start 0) (end (length-of seq)))
  (%to-string-encoded seq #'decode-rna-4bit start end :lowercase))

(defmethod coerce-sequence ((seq encoded-aa-sequence) (type (eql 'string))
                            &key (start 0) (end (length-of seq)))
  (%to-string-encoded seq #'decode-aa-7bit start end :uppercase))

(defmethod coerce-sequence ((seq virtual-dna-sequence)
                            (type (eql 'rna-sequence))
                            &key (start 0) (end (length-of seq)))
  (make-instance 'virtual-rna-sequence :length (- end start)))

(defmethod coerce-sequence ((seq virtual-rna-sequence)
                            (type (eql 'dna-sequence))
                            &key (start 0) (end (length-of seq)))
  (make-instance 'virtual-dna-sequence :length (- end start)))

(defmethod coerce-sequence ((seq encoded-dna-sequence)
                            (type (eql 'rna-sequence))
                            &key (start 0) (end (length-of seq)))
  (with-accessors ((vector vector-of))
      seq 
    (make-instance 'encoded-rna-sequence
                   :vector (token-subsequence vector start end))))

(defmethod coerce-sequence ((seq encoded-rna-sequence)
                            (type (eql 'dna-sequence))
                            &key (start 0) (end (length-of seq)))
  (with-accessors ((vector vector-of))
      seq 
    (make-instance 'encoded-dna-sequence
                   :vector (token-subsequence vector start end))))

(defmethod subsequence ((seq vector-sequence) (start fixnum)
                        &optional end)
  (with-accessors ((vector vector-of))
      seq
    (make-instance (class-of seq)
                   :vector (token-subsequence vector start end))))

(defmethod subsequence  ((seq dna-quality-sequence) (start fixnum)
                         &optional end)
  (with-accessors ((vector vector-of) (quality quality-of) (metric metric-of))
      seq
    (make-instance 'dna-quality-sequence
                   :vector (token-subsequence vector start end)
                   :quality (subseq quality start end)
                   :metric metric)))

(defmethod subsequence ((seq virtual-token-sequence) (start fixnum)
                        &optional end)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (%check-token-range length start end)
      (make-instance (class-of seq) :length (- end start)))))

(defmethod reverse-sequence ((seq vector-sequence))
  (with-accessors ((vector vector-of))
      seq 
    (make-instance (class-of seq) :vector (reverse vector))))

(defmethod reverse-sequence ((seq dna-quality-sequence))
  (with-accessors ((vector vector-of) (quality quality-of) (metric metric-of))
      seq
    (make-instance 'dna-quality-sequence
                   :vector (reverse vector)
                   :quality (reverse quality)
                   :metric metric)))

(defmethod reverse-sequence ((seq virtual-token-sequence))
  (with-accessors ((length length-of))
      seq 
    (make-instance (class-of seq) :length length)))

(defmethod nreverse-sequence ((seq vector-sequence))
  (with-accessors ((vector vector-of))
      seq
    (setf vector (nreverse vector))
    seq))

(defmethod nreverse-sequence ((seq dna-quality-sequence))
  (with-accessors ((vector vector-of) (quality quality-of))
      seq
    (setf vector (nreverse vector)
          quality (nreverse quality))
    seq))

(defmethod nreverse-sequence ((seq virtual-token-sequence))
  seq)

(defmethod complement-sequence ((seq encoded-dna-sequence))
  (with-accessors ((vector vector-of))
      seq
    (make-instance 'encoded-dna-sequence :vector
                   (complement-tokens
                    (copy-seq vector) #'complement-dna-4bit))))

(defmethod complement-sequence ((seq encoded-rna-sequence))
  (with-accessors ((vector vector-of))
      seq
    (make-instance 'encoded-rna-sequence :vector
                   (complement-tokens
                    (copy-seq vector) #'complement-rna-4bit))))

(defmethod complement-sequence ((seq virtual-token-sequence))
  (with-accessors ((length length-of))
      seq 
    (make-instance (class-of seq) :length length)))

(defmethod ncomplement-sequence ((seq encoded-dna-sequence))
  (with-accessors ((vector vector-of))
      seq
    (loop
       for i from 0 below (length vector)
       do (setf (aref vector i)
                (complement-dna-4bit (aref vector i)))))
  seq)

(defmethod ncomplement-sequence ((seq encoded-rna-sequence))
  (with-accessors ((vector vector-of))
      seq
    (loop
       for i from 0 below (length vector)
       do (setf (aref vector i)
                (complement-rna-4bit (aref vector i)))))
  seq)

(defmethod ncomplement-sequence ((seq virtual-token-sequence))
  seq)

(defmethod reverse-complement ((seq encoded-dna-sequence))
  (with-accessors ((vector vector-of))
      seq
    (make-instance 'encoded-dna-sequence :vector
                   (nreverse
                    (complement-tokens vector #'complement-dna-4bit)))))

(defmethod reverse-complement ((seq encoded-rna-sequence))
  (with-accessors ((vector vector-of))
      seq
    (make-instance 'encoded-rna-sequence :vector
                   (nreverse
                    (complement-tokens vector #'complement-rna-4bit)))))

(defmethod reverse-complement ((seq dna-quality-sequence))
  (with-accessors ((vector vector-of) (quality quality-of) (metric metric-of))
      seq
    (let ((s (make-instance 'dna-quality-sequence
                            :vector (copy-seq vector)
                            :quality (copy-seq quality)
                            :metric metric)))
      (nreverse-complement s))))

(defmethod reverse-complement ((seq virtual-token-sequence))
  (with-accessors ((length length-of))
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
  (with-accessors ((vector vector-of))
      seq
    (let ((frequencies (make-array (size-of alphabet)
                                   :element-type 'fixnum :initial-element 0)))
      (loop
         for elt across vector
         do (incf (aref frequencies (token-index alphabet elt))))
      (pairlis (copy-list (tokens-of alphabet))
               (coerce frequencies 'list)))))

(defmethod residue-position (character (seq encoded-dna-sequence)
                             &key from-end test test-not (start 0) end)
  (with-accessors ((vector vector-of))
      seq
    (%check-token-range (length vector) start end)
    (position character vector :from-end from-end :test test :test-not test-not
              :start start :end end :key #'decode-dna-4bit)))

(defmethod residue-position (character (seq encoded-rna-sequence)
                             &key from-end test test-not (start 0) end)
  (with-accessors ((vector vector-of))
      seq
    (%check-token-range (length vector) start end)
    (position character vector :from-end from-end :test test :test-not test-not
              :start start :end end :key #'decode-rna-4bit)))

(defmethod residue-position (character (seq encoded-aa-sequence)
                             &key from-end test test-not (start 0) end)
  (with-accessors ((vector vector-of))
      seq
    (%check-token-range (length vector) start end)
    (position character vector :from-end from-end :test test :test-not test-not
              :start start :end end :key #'decode-aa-7bit)))

(defmethod quality-position ((quality fixnum) (seq dna-quality-sequence)
                             &key from-end test test-not (start 0) end)
  (with-accessors ((q quality-of))
      seq
    (%check-token-range (length q) start end)
    (position quality q :from-end from-end :test test :test-not test-not
              :start start :end end)))

(defmethod translate :before ((seq na-sequence) (code genetic-code)
                              &key (start 0) end initiator-codon
                              partial-codon)
  (declare (ignore initiator-codon))
  (let* ((na-len (length-of seq))
         (end (or end na-len)))
    (%check-token-range na-len start end)
    (if (and (plusp (rem (- end start) +codon-size+))
             (not partial-codon))
        (error 'translation-error :sequence seq :start start :end end
               :genetic-code code
               :text "the translated region includes a partial codon"))))

(defmethod translate ((seq virtual-dna-sequence) (code genetic-code)
                      &key (start 0) end (initiator-codon t) partial-codon)
  (declare (ignore partial-codon))
  (translate-virtual seq code start end initiator-codon))

(defmethod translate ((seq virtual-rna-sequence) (code genetic-code)
                      &key (start 0) end (initiator-codon t) partial-codon)
  (declare (ignore partial-codon))
  (translate-virtual seq code start end initiator-codon))

(defmethod translate ((seq encoded-dna-sequence) (code genetic-code)
                      &key (start 0) end (initiator-codon t) partial-codon)
  (declare (ignore partial-codon))
  (translate-encoded-4bit seq code start end initiator-codon))

(defmethod translate ((seq encoded-rna-sequence) (code genetic-code)
                      &key (start 0) end (initiator-codon t) partial-codon)
  (declare (ignore partial-codon))
  (translate-encoded-4bit seq code start end initiator-codon))

;;; Utility functions
(defun %print-seq (name seq stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of seq)))
    (if (<= len *sequence-print-limit*)
        (format stream "#<~a \"~a\">" name (coerce-sequence seq 'string))
      (format stream "#<~a length ~d>" name len))))

(defun %print-quality-seq (name seq stream)
  "Helper function for printing bio-sequence objects."
  (with-accessors ((quality quality-of) (metric metric-of))
      seq
    (let ((len (length-of seq)))
      (if (<= len *sequence-print-limit*)
          (format stream "#<~a \"~a\" ~a quality \"~a\">"
                  name (coerce-sequence seq 'string) metric
                  (quality-string quality metric))
        (format stream "#<~a ~a quality, length ~d>"
                name metric len)))))

(defun quality-string (quality metric)
  "Wrapper for ENCODE-QUALITY that encodes QUALITY, an array of bytes
representing base quality scores, as a string using an encoder
appropriate for METRIC, a quality metric ( :PHRED or :ILLUMINA )."
  (let ((encoder (ecase metric
                   (:phred #'encode-phred-quality)
                   (:illumina #'encode-illumina-quality))))
    (encode-quality quality encoder)))

(defun encode-quality (quality encoder)
  "Encodes QUALITY, an array of bytes representing base quality
scores, as a string using function ENCODER."
  (declare (optimize (speed 3)))
  (declare (type (simple-array quality-score (*)) quality)
           (type function encoder))
  (let ((quality-str (make-string (length quality)
                                  :element-type 'base-char)))
    (copy-array quality 0 (1- (length quality))
                quality-str 0 encoder)
    quality-str))

(defun decode-quality (quality decoder)
  "Decodes the QUALITY, a string, as into a new array using function
DECODER."
  (declare (optimize (speed 3)))
  (declare (type simple-string quality)
           (type function decoder))
  (let ((quality-seq (make-array (length quality)
                                 :element-type 'quality-score)))
    (copy-array quality 0 (1- (length quality))
                quality-seq 0 decoder)
    quality-seq))

(defun ensure-encoded (vector encoder element-type)
  "If VECTOR is a string, returns an encoded token vector of element
type (unsigned-byte 4) of the same length, otherwise returns
VECTOR. ENCODER is the encoding function used to convert characters to
\(unsigned-byte 4\))."
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
    (%check-token-range (length tokens) start end)
    (let ((sub-seq (make-array (- end start)
                               :element-type
                               (array-element-type tokens))))
      (copy-array tokens start (1- end)
                  sub-seq 0)
      sub-seq)))

(defun %check-token-range (length start &optional end)
  (let ((end (or end length)))
    (cond ((or (< start 0) (> start length))
           (error 'invalid-argument-error
                  :params 'start
                  :args start
                  :text "start must be >0 and <= sequence length"))
          ((or (< end 0) (> end length))
           (error 'invalid-argument-error
                  :params 'end
                  :args end
                  :text "end must be >0 and <= sequence length"))
          ((< end start)
           (error 'invalid-argument-error
                  :params '(start end)
                  :args (list start end)
                  :text "end must be equal to or greater than start"))
          (t t))))

(defun complement-tokens (tokens comp-fn &optional (start 0) end)
  "Returns a complemented copy of TOKENS populated with elements
from TOKENS that have been transformed by COMP-FN, starting at the
first element, or index START, and continuing to the last residue, or
index END."
  (let* ((end (or end (length tokens)))
         (comp-seq (make-array (- end start)
                               :element-type
                               (array-element-type tokens))))
    (copy-array tokens start (1- end)
                comp-seq 0 comp-fn)
    comp-seq))

(defun %to-string-virtual (seq char start end token-case)
  (%check-token-range (length-of seq) start end)
  (make-string (- end start) :element-type 'base-char
               :initial-element (ecase token-case
                                  ((nil) char)
                                  (:lowercase (char-downcase char))
                                  (:uppercase char))))

(defun %to-string-encoded (seq decoder start end token-case)
  (with-accessors ((vector vector-of))
      seq
    (let ((str (make-string (- end start) :element-type 'base-char))
          (seq-end (1- end))
          (str-start 0))
      (when (< 0 (length str))
        (copy-array vector start seq-end
                    str str-start decoder))
      (ecase token-case
        ((nil) str)
        (:lowercase str)
        (:uppercase (nstring-upcase str))))))

(defmethod %ambiguousp ((seq vector-sequence) (alphabet (eql *dna*)))
  (with-accessors ((vector vector-of))
      seq
    (loop
       for elt across vector
       thereis (< 1 (length (explode-encoded-base elt))))))

(defmethod %ambiguousp ((seq vector-sequence) (alphabet (eql *rna*)))
  (with-accessors ((vector vector-of))
      seq
    (loop
       for elt across vector
       thereis (< 1 (length (explode-encoded-base elt))))))

(defmethod %ambiguousp ((seq vector-sequence) (alphabet (eql *aa*)))
  (with-accessors ((vector vector-of))
      seq
    (loop
       for elt across vector
       thereis (< 1 (length (explode-encoded-aa elt))))))

(defun %explode-ambiguity (char encoder decoder)
  "Returns a list of the ambiguity characters represented by CHAR."
  (sort (mapcar decoder (explode-encoded-base (funcall encoder char)))
        #'char<=))

(defun codon-start-p (base-pos)
  (zerop (mod base-pos +codon-size+)))

(defun codon-end-p (base-pos)
  (= 2 (mod base-pos +codon-size+)))

(defun translate-virtual (seq code start end initiator-codon)
  (let* ((na-len (length-of seq))
         (end (or end na-len))
         (aa-len (floor (- end start) +codon-size+)))
    (loop
       for i from start below end
       for j = 0 then (1+ j)
       for initiator = (and initiator-codon (< j +codon-size+))
       with codon = (mapcar #'encode-dna-4bit '(#\n #\n #\n))
       when (codon-end-p j)
       collect (translate-codon codon code
                                :initiator initiator) into aa
       finally (return (make-aa
                        (make-array aa-len
                                    :element-type '(unsigned-byte 7)
                                    :initial-contents aa))))))

(defun translate-encoded-4bit (seq code start end initiator-codon)
  (with-accessors ((vector vector-of))
      seq
    (let* ((na-len (length vector))
           (end (or end na-len))
           (aa-len (floor (- end start) +codon-size+)))
      (loop
         for i from start below end
         for j = 0 then (1+ j)
         for initiator = (and initiator-codon (< j +codon-size+))
         for codon = () then (if (codon-start-p j)
                                 ()
                               codon)
         do (push (aref vector i) codon)
         when (codon-end-p j)
         collect (translate-codon (nreverse codon) code
                                  :initiator initiator) into aa
         finally (return (make-aa
                          (make-array aa-len
                                      :element-type '(unsigned-byte 7)
                                      :initial-contents aa)))))))
