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

(defun find-alphabet (name)
  (multiple-value-bind (alphabet presentp)
      (gethash name *alphabets*)
    (unless presentp
      (error 'invalid-argument-error
             :params 'name
             :args name
             :text "no such alphabet"))
    alphabet))

(defmethod anonymousp ((seq identity-mixin))
  (null (identity-of seq)))

(defmethod size-of ((alphabet alphabet))
  (length (tokens-of alphabet)))

(defmethod token-index ((alphabet alphabet) encoded-token)
  (gethash encoded-token (index-of alphabet)))

(defmethod memberp ((alphabet alphabet) (char character))
  (contains-char-p (tokens-of alphabet) char))

(defmethod simplep ((residues string) (alphabet alphabet))
  (let ((simple (tokens-of alphabet)))
    (loop
       for residue across residues
       always (find (char-downcase residue) simple))))

(defmethod initialize-instance :after ((seq dna-sequence) &key)
  (initialize-seq seq #'ensure-encoded-4bit #'encode-dna-4bit))

(defmethod initialize-instance :after ((seq rna-sequence) &key)
  (initialize-seq seq #'ensure-encoded-4bit #'encode-rna-4bit))

(defmethod initialize-instance :after ((seq quality-mixin) &key)
  (with-slots (residues metric quality) seq
    (unless (= (length residues)
               (length quality))
      (error 'invalid-argument-error
             :params '(residues quality)
             :args (list residues quality)
             :text "the residues and quality vectors were not the same length"))
    (when (subtypep (array-element-type quality) 'character)
      (let ((decoder (ecase metric
                       (:phred #'decode-phred-quality)
                       (:illumina #'decode-illumina-quality))))
        (setf quality (decode-quality quality decoder))))))

(defmethod length-of ((seq bio-sequence))
  (with-slots (residues length) seq
    (if (null residues)
        length
      (length residues))))

(defmethod (setf length-of) (value (seq bio-sequence))
  (with-slots (residues length) seq
    (if (null residues)
        (setf length value)
      (error 'invalid-operation-error :text
             "the length of a concrete sequence may not be changed"))))

(defmethod virtualp ((seq bio-sequence))
  (null (slot-value seq 'residues)))

(defmethod residue-of :around ((seq bio-sequence) index)
  (when (virtualp seq)
    (error 'invalid-operation-error :text
           "cannot access a residue in a virtual sequence"))
  (call-next-method))

(defmethod (setf residue-of) :around (value (seq bio-sequence) index)
  (when (virtualp seq)
    (error 'invalid-operation-error :text
           "cannot access a residue in a virtual sequence"))
  (call-next-method))

(defmethod residue-of ((seq dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((residues (residues-of seq)))
    (declare (type (encoded-residues 4) residues))
    (decode-dna-4bit (aref residues index))))

(defmethod (setf residue-of) (value (seq dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((residues (residues-of seq)))
    (declare (type (encoded-residues 4) residues))
    (setf (aref residues index) (encode-dna-4bit value))))

(defmethod residue-of ((seq rna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((residues (residues-of seq)))
    (declare (type (encoded-residues 4) residues))
    (decode-rna-4bit (aref residues index))))

(defmethod (setf residue-of) (value (seq rna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((residues (residues-of seq)))
    (declare (type (encoded-residues 4) residues))
    (setf (aref residues index) (encode-rna-4bit value))))




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


(defmethod to-string :around ((seq bio-sequence)
                              &key start end token-case)
  (declare (ignore start end token-case))
  (if (virtualp seq)
      "<virtual>"
    (call-next-method)))

(defmethod to-string ((seq dna-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (let ((residues (residues-of seq))
        (str (make-string (- end start) :element-type 'base-char))
        (seq-end (1- end))
        (str-start 0))
    (declare (type (encoded-residues 4) residues)
             (type simple-base-string str)
             (type array-index seq-end str-start))
    (when (< 0 (length str))
      (copy-array residues start seq-end
                  str str-start #'decode-dna-4bit))
    (ecase token-case
      ((nil) str)
      (:lowercase str)
      (:uppercase (nstring-upcase str)))))

(defmethod to-string ((seq rna-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (let ((residues (residues-of seq))
        (str (make-string (- end start) :element-type 'base-char))
        (seq-end (1- end))
        (str-start 0))
    (declare (type (encoded-residues 4) residues)
             (type simple-base-string str)
             (type array-index seq-end str-start))
    (when (< 0 (length str))
      (copy-array residues start seq-end
                  str str-start #'decode-rna-4bit))
    (ecase token-case
      ((nil) str)
      (:lowercase str)
      (:uppercase (nstring-upcase str)))))


;;; FIXME -- add a means of making a reversed and/or complemented view
;;; of a sequence without modifying it

(defun reverse-complement-index (index len)
  (- len index))

(defmethod subsequence :around ((seq bio-sequence) (start fixnum)
                                &optional end)
  (declare (ignore start end))
  (when (virtualp seq)
    (error 'invalid-operation-error
           :text "cannot subsequence a virtual sequence"))
  (call-next-method))

(defmethod subsequence ((seq bio-sequence) (start fixnum)
                        &optional end)
  (make-instance (class-of seq)
                 :residues (residue-subsequence (residues-of seq)
                                               start end)))

(defmethod subsequence  ((seq dna-quality-sequence) (start fixnum)
                         &optional end)
  (make-instance 'dna-quality-sequence
                 :residues (residue-subsequence (residues-of seq)
                                               start end)
                 :quality (subseq (quality-of seq) start end)
                 :metric (metric-of seq)))

(defmethod reverse-sequence ((seq bio-sequence))
  (make-instance (class-of seq)
                 :residues (reverse (residues-of seq))))

(defmethod reverse-sequence ((seq dna-quality-sequence))
  (make-instance 'dna-quality-sequence
                 :residues (reverse (residues-of seq))
                 :quality (reverse (quality-of seq))
                 :metric (metric-of seq)))

(defmethod nreverse-sequence ((seq bio-sequence))
  (setf (residues-of seq) (nreverse (residues-of seq)))
  seq)

(defmethod nreverse-sequence ((seq dna-quality-sequence))
  (setf (residues-of seq) (nreverse (residues-of seq))
        (quality-of seq) (nreverse (quality-of seq)))
  seq)

(defmethod complement-sequence ((seq dna-sequence))
  (make-instance 'dna-sequence :residues
                 (complement-residues
                  (copy-seq (residues-of seq)) #'complement-dna-4bit)))

(defmethod ncomplement-sequence ((seq dna-sequence))
  (let ((residues (residues-of seq)))
    (loop
         for i from 0 below (length residues)
         do (setf (aref residues i)
                  (complement-dna-4bit (aref residues i)))))
  seq)

(defmethod reverse-complement ((seq dna-sequence))
  (make-instance 'dna-sequence :residues
                 (nreverse
                  (complement-residues
                   (residues-of seq) #'complement-dna-4bit))))

(defmethod reverse-complement ((seq dna-quality-sequence))
  (let ((s (make-instance 'dna-quality-sequence
                          :residues (copy-seq (residues-of seq))
                          :quality (copy-seq (quality-of seq))
                          :metric (metric-of seq))))
    (nreverse-complement s)))

(defmethod nreverse-complement ((seq dna-sequence))
  (nreverse-sequence (ncomplement-sequence seq)))

(defmethod nreverse-complement ((seq dna-quality-sequence))
  (nreverse-sequence (ncomplement-sequence seq)))

(defmethod search-sequence ((seq1 bio-sequence) (seq2 bio-sequence)
                            &key from-end start1 start2 end1 end2)
  (if (subtypep (class-of (alphabet-of seq1))
                (class-of (alphabet-of seq2)))
      (let ((residues1 (residues-of seq1))
            (residues2 (residues-of seq2))
            (start1 (or start1 0))
            (start2 (or start2 0)))
        (search residues1 residues2 :from-end from-end
                :start1 start1 :start2 start2
                :end1 end1 :end2 end2 :test #'eq))
    nil))

(defmethod residue-frequencies :before ((seq bio-sequence))
  (when (virtualp seq)
    (error 'invalid-operation-error
           :text (msg "cannot determine residue frequencies"
                      "of a virtual sequence"))))

(defmethod residue-frequencies ((seq bio-sequence))
  (let ((frequencies (make-array (length (tokens-of (alphabet-of seq)))
                                 :element-type 'fixnum :initial-element 0))
        (residues (residues-of seq))
        (alphabet (alphabet-of seq)))
    (loop
       for residue across residues
       do (incf (aref frequencies (token-index alphabet residue))))
    (pairlis (coerce (copy-seq (tokens-of alphabet)) 'list)
             (coerce frequencies 'list))))


(defmethod print-object ((alphabet alphabet) stream)
  (format stream "<ALPHABET ~a>" (slot-value alphabet 'name)))

(defmethod print-object ((strand sequence-strand) stream)
  (with-slots (name token number) strand
      (format stream "<SEQUENCE-STRAND ~a/~a/~a>" name token number)))

(defmethod print-object ((seq dna-sequence) stream)
  (print-seq-aux "DNA-SEQUENCE" seq stream))

(defmethod print-object ((seq rna-sequence) stream)
  (print-seq-aux "RNA-SEQUENCE" seq stream))

(defmethod print-object ((seq dna-quality-sequence) stream)
  (print-quality-seq-aux "DNA-QUALITY-SEQUENCE" seq stream))

(defun print-seq-aux (name seq stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of seq)))
    (if (and len (<= len *sequence-print-limit*)
             (not (virtualp seq)))
        (format stream "<~a \"~a\">" name (to-string seq))
      (format stream "<~a length ~d>" name len))))

(defun print-quality-seq-aux (name seq stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of seq)))
    (if (and len (<= len *sequence-print-limit*)
             (not (virtualp seq)))
        (format stream "<~a \"~a\" ~a quality \"~a\">"
                name (to-string seq) (metric-of seq)
                (quality-string (quality-of seq) (metric-of seq)))
      (format stream "<~a ~a quality, length ~d>"
              name (metric-of seq) len))))

(defun quality-string (quality metric)
  (let ((encoder (ecase metric
                   (:phred #'encode-phred-quality)
                   (:illumina #'encode-illumina-quality)))
        (str (make-string (length quality) :element-type 'base-char)))
    (map-into str encoder quality)))

(defun process-residues-args (residues length)
  "Returns its arguments, having checked their consistency for use
when making bio-sequence instances."
  (cond ((and (null residues)
              (null length))
         (error 'invalid-argument-error
                :params '(residues length)
                :args (list residues length)
                :text "expected one to be non-NIL"))
        ((and residues
              (null length))
         (unless (and (vectorp residues)
                      (not (zerop (length residues))))
           (error 'invalid-argument-error
                  :params 'residues
                  :args residues
                  :text "expected a non-empty vector"))
         (values residues (length residues)))
        ((and (null residues)
              length)
         (unless (typep length 'residues-length)
           (error 'invalid-argument-error
                  :params 'length
                  :args length
                  :text "expected a fixnum >= 1"))
         (values residues length))
        (t
         (error 'invalid-argument-error
                :params '(residues length)
                :args (list residues length)
                :text "expected one to be NIL"))))

(defun initialize-seq (seq ensure-encoded-fn encoder)
  "Initializes bio-sequence SEQ by checking using ENSURE-ENCODED-FN
that its residues are encoded."
  (with-slots (residues length) seq
    (multiple-value-bind (valid-residues valid-length)
        (process-residues-args residues length)
      (if valid-residues
          (funcall ensure-encoded-fn seq encoder)
        (setf (length-of seq) valid-length))))
  seq)

(defun ensure-encoded-4bit (seq encoder)
  "Returns SEQ after ensuring that its residues are encoded as
unsigned-byte 4 using ENCODER."
  (let ((current-residues (residues-of seq)))
    (if (equal '(unsigned-byte 4) (array-element-type current-residues))
        seq
      (let ((residues (make-array (length current-residues)
                                  :element-type '(unsigned-byte 4))))
        (declare (type (simple-array (unsigned-byte 4) *) residues)
                 (type function encoder))
        (when (< 0 (length residues))
          (copy-array current-residues 0 (1- (length residues))
                      residues 0 encoder))
        (setf (residues-of seq) residues)
        seq))))

(defun residue-subsequence (residues start end)
  "Returns a subsequence of RESIDUES between indices START and END."
  (let* ((end (or end (length residues)))
         (sub-seq (make-array (- end start)
                              :element-type
                              (array-element-type residues))))
    (copy-array residues start (1- end)
                sub-seq 0)
    sub-seq))

(defun complement-residues (residues comp-fn &optional (start 0) end)
  "Returns a complemented copy of RESIDUES populated with elements
from RESIDUES that have been transformed by COMP-FN, starting at the
first element, or index START, and continuing to the last residue, or
index END."
  (let* ((end (or end (length residues)))
         (comp-seq (make-array (- end start)
                               :element-type
                               (array-element-type residues))))
    (copy-array residues start (1- end)
                comp-seq 0 comp-fn)
    comp-seq))
