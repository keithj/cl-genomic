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

(defmethod anonymousp ((obj identity-mixin))
  (null (identity-of obj)))

(defmethod size-of ((alphabet alphabet))
  (length (tokens-of alphabet)))

(defmethod memberp ((alphabet alphabet) (char character))
  (contains-char-p (tokens-of alphabet) char))

(defmethod simplep ((token-seq string) (alphabet alphabet))
  (let ((simple (tokens-of alphabet)))
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

(defmethod slot-unbound (class (obj alphabet) (slot (eql 'decoded-index)))
  (let ((index-table (make-hash-table))
        (tokens (tokens-of obj)))
    (loop
       for i from 0 below (length tokens)
       do (setf (gethash (aref tokens i) index-table) i))
    (setf (slot-value obj 'decoded-index)
          (lambda (element)
            (gethash element index-table)))))

(defmethod initialize-instance :after ((seq dna-sequence) &key)
  (initialize-seq seq #'ensure-encoded-4bit #'encode-dna-4bit))

(defmethod initialize-instance :after ((seq rna-sequence) &key)
  (initialize-seq seq #'ensure-encoded-4bit #'encode-rna-4bit))

(defmethod initialize-instance :after ((seq quality-mixin) &key)
  (with-slots (token-seq metric quality) seq
    (unless (= (length token-seq)
               (length quality))
      (error 'invalid-argument-error
             :params '(token-seq quality)
             :args (list token-seq quality)
             :text "the token-seq and quality vectors were not the same length"))
    (when (eql 'character (array-element-type quality))
      (let ((decoder (ecase metric
                       (:phred #'decode-phred-quality)
                       (:illumina #'decode-illumina-quality))))
        (setf quality (decode-quality quality decoder))))))

(defmethod length-of ((seq bio-sequence))
  (with-slots (token-seq length) seq
    (if (null token-seq)
        length
      (length token-seq))))

(defmethod (setf length-of) (value (seq bio-sequence))
  (with-slots (token-seq length) seq
    (if (null token-seq)
        (setf length value)
      (error 'invalid-operation-error :text
             "the length of a concrete sequence may not be changed"))))

(defmethod virtualp ((seq bio-sequence))
  (null (slot-value seq 'token-seq)))

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
  (let ((token-seq (token-seq-of seq)))
    (declare (type (encoded-tokens 4) token-seq))
    (decode-dna-4bit (aref token-seq index))))

(defmethod (setf residue-of) (value (seq dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((token-seq (token-seq-of seq)))
    (declare (type (encoded-tokens 4) token-seq))
    (setf (aref token-seq index) (encode-dna-4bit value))))

(defmethod residue-of ((seq rna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((token-seq (token-seq-of seq)))
    (declare (type (encoded-tokens 4) token-seq))
    (decode-rna-4bit (aref token-seq index))))

(defmethod (setf residue-of) (value (seq rna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (let ((token-seq (token-seq-of seq)))
    (declare (type (encoded-tokens 4) token-seq))
    (setf (aref token-seq index) (encode-rna-4bit value))))




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
  (let ((token-seq (token-seq-of seq))
        (str (make-string (- end start) :element-type 'base-char))
        (seq-end (1- end))
        (str-start 0))
    (declare (type (encoded-tokens 4) token-seq)
             (type simple-base-string str)
             (type array-index seq-end str-start))
    (when (< 0 (length str))
      (copy-array token-seq start seq-end
                  str str-start #'decode-dna-4bit))
    (ecase token-case
      ((nil) str)
      (:lowercase str)
      (:uppercase (nstring-upcase str)))))

(defmethod to-string ((seq rna-sequence) &key
                      (start 0) (end (length-of seq)) token-case)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index start end))
  (let ((token-seq (token-seq-of seq))
        (str (make-string (- end start) :element-type 'base-char))
        (seq-end (1- end))
        (str-start 0))
    (declare (type (encoded-tokens 4) token-seq)
             (type simple-base-string str)
             (type array-index seq-end str-start))
    (when (< 0 (length str))
      (copy-array token-seq start seq-end
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
                 :token-seq (token-subsequence (token-seq-of seq)
                                               start end)))

(defmethod subsequence  ((seq dna-quality-sequence) (start fixnum)
                         &optional end)
  (make-instance 'dna-quality-sequence
                 :token-seq (token-subsequence (token-seq-of seq)
                                               start end)
                 :quality (subseq (quality-of seq) start end)
                 :metric (metric-of seq)))

(defmethod reverse-sequence ((seq bio-sequence))
  (make-instance (class-of seq)
                 :token-seq (reverse (token-seq-of seq))))

(defmethod reverse-sequence ((seq dna-quality-sequence))
  (make-instance 'dna-quality-sequence
                 :token-seq (reverse (token-seq-of seq))
                 :quality (reverse (quality-of seq))
                 :metric (metric-of seq)))

(defmethod nreverse-sequence ((seq bio-sequence))
  (setf (token-seq-of seq) (nreverse (token-seq-of seq)))
  seq)

(defmethod nreverse-sequence ((seq dna-quality-sequence))
  (setf (token-seq-of seq) (nreverse (token-seq-of seq))
        (quality-of seq) (nreverse (quality-of seq)))
  seq)

(defmethod complement-sequence ((seq dna-sequence))
  (make-instance 'dna-sequence :token-seq
                 (complement-token-seq
                  (copy-seq (token-seq-of seq)) #'complement-dna-4bit)))

(defmethod ncomplement-sequence ((seq dna-sequence))
  (let ((token-seq (token-seq-of seq)))
    (loop
         for i from 0 below (length token-seq)
         do (setf (aref token-seq i)
                  (complement-dna-4bit (aref token-seq i)))))
  seq)

(defmethod reverse-complement ((seq dna-sequence))
  (make-instance 'dna-sequence :token-seq
                 (nreverse
                  (complement-token-seq
                   (token-seq-of seq) #'complement-dna-4bit))))

(defmethod reverse-complement ((seq dna-quality-sequence))
  (let ((s (make-instance 'dna-quality-sequence
                          :token-seq (copy-seq (token-seq-of seq))
                          :quality (copy-seq (quality-of seq))
                          :metric (metric-of seq))))
    (nreverse-complement s)))

(defmethod nreverse-complement ((seq dna-sequence))
  (nreverse-sequence (ncomplement-sequence seq)))

(defmethod nreverse-complement ((seq dna-quality-sequence))
  (nreverse-sequence (ncomplement-sequence seq)))

(defmethod residue-frequencies :before ((seq bio-sequence))
  (when (virtualp seq)
    (error 'invalid-operation-error
           :text (msg "cannot determine residue frequencies"
                      "of a virtual sequence"))))

(defmethod search-sequence ((seq1 bio-sequence) (seq2 bio-sequence)
                            &key from-end start1 start2 end1 end2)
  (if (subtypep (class-of (alphabet-of seq1))
                (class-of (alphabet-of seq2)))
      (let ((token-seq1 (token-seq-of seq1))
            (token-seq2 (token-seq-of seq2))
            (start1 (or start1 0))
            (start2 (or start2 0)))
        (search token-seq1 token-seq2 :from-end from-end
                :start1 start1 :start2 start2 :end1 end1 :end2 end2 :test #'eq))
    nil))



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

(defmethod print-object ((obj dna-sequence) stream)
  (print-seq-aux "DNA-SEQUENCE" obj stream))

(defmethod print-object ((obj rna-sequence) stream)
  (print-seq-aux "RNA-SEQUENCE" obj stream))

(defmethod print-object ((obj dna-quality-sequence) stream)
  (print-quality-seq-aux "DNA-QUALITY-SEQUENCE" obj stream))

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
        (format stream "<~a \"~a\" ~a quality \"~a\">"
                name (to-string obj) (metric-of obj)
                (quality-string (quality-of obj) (metric-of obj)))
      (format stream "<~a ~a quality, length ~d>"
              name (metric-of obj) len))))

(defun quality-string (quality metric)
  (let ((encoder (ecase metric
                   (:phred #'encode-phred-quality)
                   (:illumina #'encode-illumina-quality)))
        (str (make-string (length quality) :element-type 'base-char)))
    (map-into str encoder quality)))

(defun process-token-seq-args (token-seq length)
  "Returns its arguments, having checked their consistency for use
when making bio-sequence instances."
  (cond ((and (null token-seq)
              (null length))
         (error 'invalid-argument-error
                :params '(token-seq length)
                :args (list token-seq length)
                :text "expected one to be non-NIL"))
        ((and token-seq
              (null length))
         (unless (and (vectorp token-seq)
                      (not (zerop (length token-seq))))
           (error 'invalid-argument-error
                  :params 'token-seq
                  :args token-seq
                  :text "expected a non-empty vector"))
         (values token-seq (length token-seq)))
        ((and (null token-seq)
              length)
         (unless (typep length 'token-seq-length)
           (error 'invalid-argument-error
                  :params 'length
                  :args length
                  :text "expected a fixnum >= 1"))
         (values token-seq length))
        (t
         (error 'invalid-argument-error
                :params '(token-seq length)
                :args (list token-seq length)
                :text "expected one to be NIL"))))

(defun initialize-seq (bio-seq ensure-encoded-fn encoder)
  "Initializes bio-sequence BIO-SEQ by checking using
ENSURE-ENCODED-FN that its tokens are encoded."
  (with-slots (token-seq length) bio-seq
    (multiple-value-bind (valid-token-seq valid-length)
        (process-token-seq-args token-seq length)
      (if valid-token-seq
          (funcall ensure-encoded-fn bio-seq encoder)
        (setf (length-of bio-seq) valid-length))))
  bio-seq)

(defun ensure-encoded-4bit (bio-seq encoder)
  "Returns BIO-SEQ after ensuring that its tokens are encoded as
unsigned-byte 4 using ENCODER."
  (let ((current-token-seq (token-seq-of bio-seq)))
    (if (equal '(unsigned-byte 4) (array-element-type current-token-seq))
        bio-seq
      (let ((token-seq (make-array (length current-token-seq)
                                   :element-type '(unsigned-byte 4))))
        (declare (type (simple-array (unsigned-byte 4) *) token-seq)
                 (type function encoder))
        (when (< 0 (length token-seq))
          (copy-array current-token-seq 0 (1- (length token-seq))
                      token-seq 0 encoder))
        (setf (token-seq-of bio-seq) token-seq)
        bio-seq))))

(defun token-subsequence (token-seq start end)
  "Returns a subsequence of TOKEN-SEQ between indices START and END."
  (let* ((end (or end (length token-seq)))
         (sub-seq (make-array (- end start)
                              :element-type
                              (array-element-type token-seq))))
    (copy-array token-seq start (1- end)
                sub-seq 0)
    sub-seq))

(defun complement-token-seq (token-seq comp-fn &optional (start 0) end)
  "Returns a complemented copy of TOKEN-SEQ populated with elements
from TOKEN-SEQ that have been transformed by COMP-FN, starting at the
first element, or index START, and continuing to the last residue, or
index END."
  (let* ((end (or end (length token-seq)))
         (comp-seq (make-array (- end start)
                               :element-type
                               (array-element-type token-seq))))
    (copy-array token-seq start (1- end)
                comp-seq 0 comp-fn)
    comp-seq))
