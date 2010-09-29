;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
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
  "Parameterized type for encoded sequence residues."
  `(simple-array (unsigned-byte ,n) (*)))

(deftype quality-score ()
  "Type for sequence base quality score."
  '(signed-byte 8))

(defparameter *sequence-print-limit* 160
  "Maximum length of sequence to be pretty-printed.")

(defmacro define-strand-decoder (type test forward reverse unknown)
  "Defines a sequence strand decoder method for a particular strand
representation that converts an ad hoc strand designator to a
canonical strand object. There is no standard nomenclature for
sequence strands and they may be represented as strings, characters,
integers, mathematical symbols etc. For example +/-/? or fwd/rev/unk
or 1/-1/0.

Arguments:

- type (Lisp type): The Lisp type of the strand designators.
- test (function): The test function to be used to compare strand
  designators.
- forward (object): The forward strand designator.
- reverse (object): The reverse strand designator.
- unknown (object): The unknown strand designator."
  `(progn
    (defmethod decode-strand ((strand ,type) &key strict)
      (cond ((,test ,forward strand)
             *forward-strand*)
            ((,test ,reverse strand)
             *reverse-strand*)
            ((,test ,unknown strand)
             *unknown-strand*)
            (t
             (when strict
               (check-arguments nil (strand)
                                "not a valid ~a strand designator" ',type)))))))

(defmacro with-sequence-residues ((var seq &key start end) &body body)
  "Executes BODY in the context of bio-sequence SEQ such that VAR is
iteratively bound to each residue of SEQ in turn.

Arguments:

- var (symbol or symbols): The iteration variable to which the current
residue will be bound during iteration.
- seq (bio-sequence): The sequence iterated to be iterated over.

Key:

- start (fixnum): The start index of the iteration.
- end (fixnum): The end index of the iteration.

Body:

Forms to be executed in the context of these bindings."
  (with-gensyms (lower upper vector)
    `(let ((,lower (or ,start 0))
           (,upper (or ,end (length-of ,seq))))
      (etypecase ,seq
        (encoded-dna-sequence
         (with-slots ((,vector vector))
             ,seq
           (loop
              for i of-type vector-index from ,lower below ,upper
              do (symbol-macrolet ((,var (%aref-dna-4bit ,vector i)))
                   ,@body))))
        (t
         (loop
            for i of-type vector-index from ,lower below ,upper
            do (symbol-macrolet ((,var (element-of ,seq i)))
                 ,@body))))
      ,seq)))

(define-strand-decoder string string= "+" "-" "?")
(define-strand-decoder fixnum = 1 -1 0)
(define-strand-decoder character char= #\+ #\- #\?)
(define-strand-decoder symbol eql :forward :reverse :unknown)

(defun make-dna (residues &rest initargs)
  "Returns a new DNA sequence object.

Arguments:

- residues (vector or NIL): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.

Rest:

- initargs: Any initialization arguments.

Returns:

- A DNA sequence object."
  (if (or (null residues) (zerop (length residues)))
      (apply #'make-instance 'virtual-dna-sequence initargs)
      (make-encoded-vector-seq 'encoded-dna-sequence residues
                               #'encode-dna-4bit '(unsigned-byte 4) initargs)))

(defun make-rna (residues &rest initargs)
  "Returns a new RNA sequence object.

Arguments:

- residues (vector or NIL): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.

Rest:

- initargs: Any initialization arguments.

Returns:

- A RNA sequence object."
  (if (or (null residues) (zerop (length residues)))
      (apply #'make-instance 'virtual-rna-sequence initargs)
      (make-encoded-vector-seq 'encoded-rna-sequence residues
                               #'encode-rna-4bit '(unsigned-byte 4) initargs)))

(defun make-dna-quality (residues quality &rest initargs
                         &key (metric :sanger) &allow-other-keys)
  "Returns a new DNA sequence object with quality.

Arguments:

- residues (vector): A vector of characters or encoded residues. It is
  not possible to create virtual sequences with quality information.
- quality (vector integer): A vector of integers, the same length as
  RESIDUES, containing quality data.

Rest:

- initargs: Any initialization arguments.

Key:

- metric (symbol): A symbol indicating the quality metric used in
  QUALITY. Either :sanger , :solexa or :illumina.

Returns:

- A DNA sequence object."
  (check-arguments (= (length residues) (length quality)) (residues quality)
                   "the residues and quality vectors were not the same length")
  (let ((initargs (remove-key-values '(:metric) initargs)))
    (apply #'make-instance 'dna-quality-sequence
           :vector (ensure-encoded
                    residues #'encode-dna-4bit '(unsigned-byte 4))
           :quality (ensure-decoded-quality quality metric)
           :metric metric
           initargs)))

(defun make-aa (residues &rest initargs)
  "Returns a new AA sequence object.

Arguments:

- residues (vector): A vector of characters, encoded residues
  or NIL. The latter creates a virtual sequence and requires a length
  argument to be supplied.

Rest:

- initargs: Any initialization arguments.

Returns:

- An AA sequence object."
  (if (null residues)
      (apply #'make-instance 'virtual-aa-sequence initargs)
    (make-encoded-vector-seq 'encoded-aa-sequence residues #'encode-aa-7bit
                             '(unsigned-byte 7) initargs)))

(defun make-encoded-vector-seq (class residues encoder element-type initargs)
  (check-arguments (vectorp residues) (residues) "expected a vector")
  (apply #'make-instance class
         :vector (ensure-encoded residues encoder element-type) initargs))

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
  (check-arguments (apply #'same-biotype-p seqs) (seqs)
                   "expected all sequences to be one of DNA, RNA or AA")
  (let ((construct (cond ((dna-sequence-p (first seqs))
                          #'make-dna)
                         ((rna-sequence-p (first seqs))
                          #'make-rna)
                         ((aa-sequence-p (first seqs))
                          #'make-aa)
                         (t
                          (check-arguments nil (seqs)
                                           (txt "expected all sequences to be"
                                                "one of DNA, RNA or AA"))))))
    ;; We use coerce-sequence to avoid a special case for virtual
    ;; sequences. We can't simply concatenate the encoded residue
    ;; vectors because virtual sequences do not have them
    (funcall construct (apply #'concatenate 'string
                              (mapcar (lambda (s)
                                        (coerce-sequence s 'string))
                                      seqs)))))

;;; Initialization methods
(defmethod initialize-instance :after ((seq na-sequence) &key)
  (with-accessors ((num-strands num-strands-of))
      seq
    (check-arguments (< 0 num-strands 3) (num-strands)
                     "nucleic acid sequences may have 1 or 2 strands")))

;;; Printing methods
(defmethod print-object ((alphabet alphabet) stream)
  (print-unreadable-object (alphabet stream :type t)
    (princ (slot-value alphabet 'name) stream)))

(defmethod print-object ((strand sequence-strand) stream)
  (print-unreadable-object (strand stream :type t)
    (with-slots (name token number)
        strand
      (format stream "~a/~a/~a" name token number))))

(defmethod print-object ((seq dna-sequence) stream)
  (%print-seq seq stream))

(defmethod print-object ((seq dna-sequence) stream)
  (%print-seq seq stream))

(defmethod print-object ((seq rna-sequence) stream)
  (%print-seq seq stream))

(defmethod print-object ((seq dna-quality-sequence) stream)
  (%print-quality-seq seq stream))

(defmethod print-object ((seq aa-sequence) stream)
  (%print-seq seq stream))

(defmethod print-object ((seq mapped-dna-sequence) stream)
  (if (dxn:in-memory-p seq)
      (%print-seq seq stream)
    (print-unreadable-object (seq stream :type t :identity t)
      (format stream "length ~d UNMAPPED" (length-of seq)))))

;;; Implementation methods
(defmethod anonymousp ((seq identity-mixin))
  (null (identity-of seq)))

(defmethod strand-designator-p (strand)
  (declare (ignore strand))
  nil)

(defmethod strand-designator-p ((strand string))
  (decode-strand strand))

(defmethod strand-designator-p ((strand fixnum))
  (decode-strand strand))

(defmethod strand-designator-p ((strand character))
  (decode-strand strand))

(defmethod strand-designator-p ((strand symbol))
  (decode-strand strand))

(defmethod forward-strand-p ((strand sequence-strand))
  (if (eql *forward-strand* strand)
      t
    nil))

(defmethod reverse-strand-p ((strand sequence-strand))
  (if (eql *reverse-strand* strand)
      t
    nil))

(defmethod unknown-strand-p ((strand sequence-strand))
  (if (eql *unknown-strand* strand)
      t
    nil))

(defmethod strand= ((strand1 sequence-strand) (strand2 sequence-strand))
  (cond ((or (eql *unknown-strand* strand1) (eql *unknown-strand* strand2))
         nil)
        (t
         (eql strand1 strand2))))

(defmethod strand-equal ((strand1 sequence-strand) (strand2 sequence-strand))
  (cond ((or (eql *unknown-strand* strand1) (eql *unknown-strand* strand2))
         t)
        (t
         (eql strand1 strand2))))

(defmethod complement-strand ((strand sequence-strand))
  (cond ((eql *forward-strand* strand)
         *reverse-strand*)
        ((eql *reverse-strand* strand)
         *forward-strand*)
        (t
         *unknown-strand*)))

(defmethod complementp ((strand1 sequence-strand) (strand2 sequence-strand))
  (strand= (complement-strand strand1) strand2))

(defmethod ambiguousp ((seq virtual-token-sequence))
  t)

;;; FIXME -- this is only implemented for DNA/RNA/AA. If a user
;;; extends the set of alphabets this must be changed.
(defmethod ambiguousp ((seq encoded-vector-sequence))
  (with-slots (vector alphabet)
      seq
    (let ((fn (cond ((eql *dna* alphabet)
                     #'enum-encoded-base)
                    ((eql *rna* alphabet)
                     #'enum-encoded-base)
                    ((eql *aa* alphabet)
                     #'enum-encoded-aa)
                    (t
                     (error "Unknown alphabet ~a" alphabet)))))
      (loop
         for elt across vector
         thereis (< 1 (length (funcall fn elt)))))))

(defmethod ambiguousp ((seq mapped-dna-sequence))
  (with-slots ((area dxn:mmap-area) length)
      seq
    (loop
       for i from 0 below length
       thereis (< 1 (length
                     (enum-dna-base
                      (code-char
                       (cffi:mem-aref (dxn:mmap-area-ptr area) :char i))))))))

(defmethod simplep ((seq virtual-token-sequence))
  nil)

(defmethod simplep ((seq encoded-vector-sequence))
  (not (ambiguousp seq)))

(defmethod simplep ((seq mapped-dna-sequence))
  (not (ambiguousp seq)))

(defmethod virtualp ((seq token-sequence))
  (declare (ignore seq))
  nil)

(defmethod virtualp ((seq virtual-token-sequence))
  (declare (ignore seq))
  t)

(defmethod length-of ((seq encoded-vector-sequence))
  (length (slot-value seq 'vector)))

(defmethod length-of ((seq mapped-dna-sequence))
  (dxn:length-of seq))

(defmethod single-stranded-p ((seq na-sequence))
  (with-slots (num-strands)
      seq
    (= 1 num-strands)))

(defmethod double-stranded-p ((seq na-sequence))
  (with-slots (num-strands)
      seq
    (= 2 num-strands)))

(defmethod (setf num-strands-of) :before (value (seq na-sequence))
  (check-arguments (< 0 value 3) (value)
                   "nucleic acid sequences may have 1 or 2 strands"))

(defmethod num-strands-of ((seq aa-sequence))
  (error 'bio-sequence-op-error
         :text "Amino-acid sequences are not stranded."))

(defmethod element-of ((seq encoded-dna-sequence) (index fixnum))
  (%aref-dna-4bit (slot-value seq 'vector) index))

(defmethod (setf element-of) (value (seq encoded-dna-sequence) (index fixnum))
  (setf (%aref-dna-4bit (slot-value seq 'vector) index) value))

(defmethod element-of ((seq encoded-rna-sequence) (index fixnum))
  (%aref-rna-4bit (slot-value seq 'vector) index))

(defmethod (setf element-of) (value (seq encoded-rna-sequence) (index fixnum))
  (setf (%aref-rna-4bit (slot-value seq 'vector) index) value))

(defmethod element-of :before ((seq virtual-token-sequence) (index fixnum))
  (%check-token-range (1- (length-of seq)) index index))

(defmethod element-of ((seq virtual-dna-sequence) (index fixnum))
  (%check-token-range (1- (length-of seq)) index index)
  #\n)

(defmethod element-of ((seq virtual-rna-sequence) (index fixnum))
  (%check-token-range (1- (length-of seq)) index index)
  #\n)

(defmethod element-of ((seq virtual-aa-sequence) (index fixnum))
  (%check-token-range (1- (length-of seq)) index index)
  #\X)

(defmethod element-of ((seq mapped-dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (with-slots ((area dxn:mmap-area))
      seq
    (code-char (cffi:mem-aref (dxn:mmap-area-ptr area) :char index))))

(defmethod (setf element-of) (value (seq mapped-dna-sequence) (index fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (with-slots ((area dxn:mmap-area))
      seq
    (declare (type fixnum index))
    (setf (cffi:mem-aref (dxn:mmap-area-ptr area) :char index)
          (char-code value))))

(declaim (inline residue-of))
(defun residue-of (seq index)
  "Returns the residue of SEQ at INDEX. This is a synonym of ELEMENT-OF."
  (element-of seq index))

(defun (setf residue-of) (value seq index)
  "Sets the residue of SEQ at INDEX to VALUE. This is a synonym of
(SETF ELEMENT-OF)."
  (setf (element-of seq index) value))

(defmethod num-gaps-of ((seq encoded-vector-sequence) &key (start 0) end)
  (let ((end (or end (length-of seq))))
    (loop
       for i from start below end
       count (= +encoded-gap-char+ (aref (slot-value seq 'vector) i)))))

(defmethod coerce-sequence ((seq dna-sequence) (type (eql 'dna-sequence))
                            &key (start 0) end)
  (maybe-subsequence seq start end))

(defmethod coerce-sequence ((seq rna-sequence) (type (eql 'rna-sequence))
                            &key (start 0) end)
  (maybe-subsequence seq start end))

(defmethod coerce-sequence ((seq aa-sequence) (type (eql 'aa-sequence))
                            &key (start 0) end)
  (maybe-subsequence seq start end))

(defmethod coerce-sequence ((seq virtual-dna-sequence) (type (eql 'string))
                            &key (start 0) end)
  (%to-string-virtual seq #\n start end :lower))

(defmethod coerce-sequence ((seq virtual-rna-sequence) (type (eql 'string))
                            &key (start 0) end)
  (%to-string-virtual seq #\n start end :lower))

(defmethod coerce-sequence ((seq virtual-aa-sequence) (type (eql 'string))
                            &key (start 0) end)
  (%to-string-virtual seq #\X start end :upper))

(defmethod coerce-sequence ((seq encoded-dna-sequence) (type (eql 'string))
                            &key (start 0) end)
  (%to-string-encoded seq #'decode-dna-4bit start end :lower))

(defmethod coerce-sequence ((seq encoded-rna-sequence) (type (eql 'string))
                            &key (start 0) end)
  (%to-string-encoded seq #'decode-rna-4bit start end :lower))

(defmethod coerce-sequence ((seq encoded-aa-sequence) (type (eql 'string))
                            &key (start 0) end)
  (%to-string-encoded seq #'decode-aa-7bit start end :upper))

(defmethod coerce-sequence ((seq virtual-dna-sequence)
                            (type (eql 'rna-sequence))
                            &key (start 0) end)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (%check-token-range length start end)
      (make-instance 'virtual-rna-sequence :length (- end start)))))

(defmethod coerce-sequence ((seq virtual-rna-sequence)
                            (type (eql 'dna-sequence))
                            &key (start 0) end)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (%check-token-range length start end)
      (make-instance 'virtual-dna-sequence :length (- end start)))))

(defmethod coerce-sequence ((seq encoded-dna-sequence)
                            (type (eql 'rna-sequence))
                            &key (start 0) end)
  (with-slots (vector)
      seq
    (let ((end (or end (length vector))))
      (make-instance 'encoded-rna-sequence
                     :vector (token-subsequence vector start end)))))

(defmethod coerce-sequence ((seq encoded-rna-sequence)
                            (type (eql 'dna-sequence))
                            &key (start 0) end)
  (with-slots (vector)
      seq
     (let ((end (or end (length vector))))
       (make-instance 'encoded-dna-sequence
                      :vector (token-subsequence vector start end)))))

(defmethod coerce-sequence ((seq mapped-dna-sequence) (type (eql 'string))
                            &key (start 0) end)
  (with-slots ((length dxn:length) (area dxn:mmap-area))
      seq
    (let* ((end (or end length))
           (str (make-array (- end start) :element-type 'base-char
                            :adjustable t :fill-pointer 0))
           (ptr (dxn:mmap-area-ptr area)))
      (if (dxn:mmap-area-live-p area) ; Only access residues if mapped
          (with-output-to-string (stream str :element-type 'base-char)
            (loop
               for i from start below end
               do (write-char (code-char (cffi:mem-aref ptr :char i)) stream)
               finally (return str)))
        (error 'invalid-operation-error
               :format-control (txt "cannot coerce ~a to a"
                                    "string when unmapped from memory")
               :format-arguments (list seq))))))

(defmethod subsequence ((seq encoded-vector-sequence) (start fixnum)
                        &optional end)
  (with-slots (vector)
      seq
    (let ((end (or end (length vector))))
      (make-instance (class-of seq)
                     :vector (token-subsequence vector start end)))))

(defmethod subsequence  ((seq dna-quality-sequence) (start fixnum)
                         &optional end)
  (with-slots (vector quality metric)
      seq
    (let ((end (or end (length vector))))
      (make-instance 'dna-quality-sequence
                     :vector (token-subsequence vector start end)
                     :quality (subseq quality start end) :metric metric))))

(defmethod subsequence ((seq virtual-token-sequence) (start fixnum)
                        &optional end)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (%check-token-range length start end)
      (make-instance (class-of seq) :length (- end start)))))

(defmethod subsequence ((seq mapped-dna-sequence) (start fixnum)
                        &optional end)
  (declare (optimize (speed 3)))
  (with-slots ((length dxn:length) (area dxn:mmap-area) )
      seq
    (let ((end (or end length)))
      (declare (fixnum end))
      (let ((vector (make-array (- end start) :element-type '(unsigned-byte 4)))
            (ptr (dxn:mmap-area-ptr area)))
        (loop
           for i of-type vector-index from start below end
           for j of-type vector-index = 0 then (1+ j)
           do (setf (aref vector j)
                    (encode-dna-4bit (code-char (cffi:mem-aref ptr :char i)))))
        (make-instance 'encoded-dna-sequence :vector vector)))))

(defmethod reverse-sequence ((seq encoded-vector-sequence))
  (make-instance (class-of seq) :vector (reverse (slot-value seq 'vector))))

(defmethod reverse-sequence ((seq dna-quality-sequence))
  (with-slots (vector quality metric)
      seq
    (make-instance 'dna-quality-sequence :vector (reverse vector)
                   :quality (reverse quality) :metric metric)))

(defmethod reverse-sequence ((seq virtual-token-sequence))
  (make-instance (class-of seq) :length (length-of seq)))

(defmethod reverse-sequence ((seq mapped-dna-sequence))
  (declare (optimize (speed 3) (safety 1)))
  (with-slots ((length dxn:length) (area dxn:mmap-area))
      seq
    (let ((rev-seq (make-instance 'mapped-dna-sequence :length length))
          (ptr (dxn:mmap-area-ptr area)))
      (with-slots ((rev-area dxn:mmap-area))
          rev-seq
        (let ((end (1- (the fixnum length)))
              (rev-ptr (dxn:mmap-area-ptr rev-area)))
          (loop
             for i from 0 to end
             do (let ((j (- end i)))
                  (setf (cffi:mem-aref rev-ptr :char i)
                        (cffi:mem-aref ptr :char j))))))
      rev-seq)))

(defmethod nreverse-sequence ((seq encoded-vector-sequence))
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

(defmethod nreverse-sequence ((seq mapped-dna-sequence))
  (declare (optimize (speed 3) (safety 1)))
  (with-slots ((length dxn:length) (area dxn:mmap-area))
      seq
    (let ((rev-seq (make-instance 'mapped-dna-sequence :length length))
          (ptr (dxn:mmap-area-ptr area)))
      (with-slots ((rev-area dxn:mmap-area))
          rev-seq
        (let ((end (1- (the fixnum length)))
              (rev-ptr (dxn:mmap-area-ptr rev-area)))
          (loop
             for i from 0 to end
             do (let ((j (- end i)))
                  (setf (cffi:mem-aref rev-ptr :char i)
                        (cffi:mem-aref ptr :char j)))))
        (dxn:free-mapped-vector seq)    ; Free the original seq
        (setf area rev-area))           ; Swap in the reversed seq
      seq)))

(defmethod complement-sequence ((seq encoded-dna-sequence))
  (make-instance 'encoded-dna-sequence
                 :vector (complement-tokens (copy-seq (slot-value seq 'vector))
                                            #'complement-dna-4bit)))

(defmethod complement-sequence ((seq encoded-rna-sequence))
  (make-instance 'encoded-rna-sequence
                 :vector (complement-tokens (copy-seq (slot-value seq 'vector))
                                            #'complement-rna-4bit)))

(defmethod complement-sequence ((seq virtual-token-sequence))
  (make-instance (class-of seq) :length (length-of seq)))

(defmethod complement-sequence ((seq mapped-dna-sequence))
  (declare (optimize (speed 3) (safety 1)))
  (with-slots ((length dxn:length) (area dxn:mmap-area))
      seq
    (let ((cmp-seq (make-instance 'mapped-dna-sequence :length length))
          (ptr (dxn:mmap-area-ptr area)))
      (with-slots ((cmp-area dxn:mmap-area))
          cmp-seq
        (let ((end (1- (the fixnum length)))
              (cmp-ptr (dxn:mmap-area-ptr cmp-area)))
          (loop
             for i from 0 to end
             do (setf (cffi:mem-aref cmp-ptr :char i)
                      (complement-dna-8bit (cffi:mem-aref ptr :char i))))))
      cmp-seq)))

(defmethod ncomplement-sequence ((seq encoded-dna-sequence))
  (with-slots (vector)
      seq
    (loop
       for i from 0 below (length vector)
       do (setf (aref vector i) (complement-dna-4bit (aref vector i)))))
  seq)

(defmethod ncomplement-sequence ((seq encoded-rna-sequence))
  (with-slots (vector)
      seq
    (loop
       for i from 0 below (length vector)
       do (setf (aref vector i) (complement-rna-4bit (aref vector i)))))
  seq)

(defmethod ncomplement-sequence ((seq virtual-token-sequence))
  seq)

(defmethod ncomplement-sequence ((seq mapped-dna-sequence))
  (declare (optimize (speed 3) (safety 1)))
  (with-slots ((length dxn:length) (area dxn:mmap-area))
      seq
    (let ((ptr (dxn:mmap-area-ptr area))
          (end (1- (the fixnum length))))
      (loop
         for i from 0 to end
         do (setf (cffi:mem-aref ptr :char i)
                  (complement-dna-8bit (cffi:mem-aref ptr :char i))))
      seq)))

(defmethod reverse-complement ((seq encoded-dna-sequence))
  (make-instance 'encoded-dna-sequence
                 :vector (nreverse (complement-tokens (slot-value seq 'vector)
                                                      #'complement-dna-4bit))))

(defmethod reverse-complement ((seq encoded-rna-sequence))
  (make-instance 'encoded-rna-sequence
                 :vector (nreverse (complement-tokens (slot-value seq 'vector)
                                                      #'complement-rna-4bit))))

(defmethod reverse-complement ((seq dna-quality-sequence))
  (with-slots (vector quality metric)
      seq
    (let ((s (make-instance 'dna-quality-sequence :vector (copy-seq vector)
                            :quality (copy-seq quality) :metric metric)))
      (nreverse-complement s))))

(defmethod reverse-complement ((seq virtual-token-sequence))
  (make-instance (class-of seq) :length (length-of seq)))

(defmethod reverse-complement ((seq mapped-dna-sequence))
  (ncomplement-sequence (reverse-sequence seq)))

(defmethod nreverse-complement ((seq token-sequence))
  (nreverse-sequence (ncomplement-sequence seq)))

(defmethod nreverse-complement ((seq virtual-token-sequence))
  seq)

(defmethod search-sequence :around ((seq1 token-sequence)
                                    (seq2 token-sequence)
                                    &key from-end start1 end1 start2 end2)
  (declare (ignore from-end start1 start2 end1 end2))
  (if (subtypep (class-of (alphabet-of seq1))
                (class-of (alphabet-of seq2)))
      (call-next-method)
    nil))

(defmethod search-sequence ((seq1 encoded-vector-sequence)
                            (seq2 encoded-vector-sequence)
                            &key from-end (start1 0) end1 (start2 0) end2)
  (let ((vector1 (slot-value seq1 'vector))
        (vector2 (slot-value seq2 'vector)))
    (search vector1 vector2 :from-end from-end :start1 start1 :start2 start2
            :end1 end1 :end2 end2 :test #'eq)))

(defmethod residue-frequencies ((seq encoded-vector-sequence))
  (with-slots (alphabet vector)
      seq
    (let ((frequencies (make-array (size-of alphabet)
                                   :element-type 'fixnum :initial-element 0)))
      (loop
         for token across vector
         do (incf (aref frequencies (token-index token alphabet))))
      (pairlis (copy-list (tokens-of alphabet)) (coerce frequencies 'list)))))

(defmethod residue-frequencies ((seq virtual-dna-sequence))
  ;; Returns 'n' x the sequence length
  (pairlis (list #\n) (list (length-of seq))))

(defmethod residue-frequencies ((seq virtual-rna-sequence))
  ;; Returns 'n' x the sequence length
  (pairlis (list #\n) (list (length-of seq))))

(defmethod residue-frequencies ((seq virtual-aa-sequence))
  ;; Returns 'X' x the sequence length
  (pairlis (list #\X) (list (length-of seq))))

(defmethod residue-frequencies ((seq mapped-dna-sequence))
  (with-slots (alphabet (length dxn:length) (area dxn:mmap-area))
      seq
    (let ((frequencies (make-array (size-of alphabet)
                                   :element-type 'fixnum :initial-element 0)))
      (loop
         for i of-type vector-index from 0 below length
         do (let ((token (encode-dna-4bit
                          (code-char
                           (cffi:mem-aref (dxn:mmap-area-ptr area)
                                          :char i)))))
              (incf (aref frequencies (token-index token alphabet))))))))

(defmethod residue-position (character (seq encoded-dna-sequence)
                             &key from-end test test-not (start 0) end)
  (with-slots (vector)
      seq
    (%check-token-range (length vector) start end)
    (position character vector
              :from-end from-end :test test :test-not test-not
              :start start :end end :key #'decode-dna-4bit)))

(defmethod residue-position (character (seq encoded-rna-sequence)
                             &key from-end test test-not (start 0) end)
  (with-slots (vector)
      seq
    (%check-token-range (length vector) start end)
    (position character vector :from-end from-end :test test :test-not test-not
              :start start :end end :key #'decode-rna-4bit)))

(defmethod residue-position (character (seq encoded-aa-sequence)
                             &key from-end test test-not (start 0) end)
  (with-slots (vector)
      seq
    (%check-token-range (length vector) start end)
    (position character vector :from-end from-end :test test :test-not test-not
              :start start :end end :key #'decode-aa-7bit)))

(defmethod quality-position ((quality fixnum) (seq dna-quality-sequence)
                             &key from-end test test-not (start 0) end)
  (with-slots ((qvector quality))
      seq
    (%check-token-range (length qvector) start end)
    (position quality qvector :from-end from-end :test test :test-not test-not
              :start start :end end)))

(defmethod translate :before ((seq na-sequence) (code genetic-code)
                              &key (start 0) end initiator-codon partial-codon)
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
(defun %print-seq (seq stream)
  "Helper function for printing bio-sequence objects."
  (print-unreadable-object (seq stream :type t :identity t)
    (let ((len (length-of seq)))
      (if (<= len *sequence-print-limit*)
          (format stream "~s" (coerce-sequence seq 'string))
        (format stream "length ~d" len)))))

(defun %print-quality-seq (seq stream)
  "Helper function for printing bio-sequence objects."
  (print-unreadable-object (seq stream :type t :identity t)
    (with-accessors ((length length-of))
        seq
      (with-slots (quality metric)
          seq
        (if (<= length *sequence-print-limit*)
            (format stream "~s ~a quality \"~a\"" (coerce-sequence seq 'string)
                    metric (quality-string quality metric))
          (format stream "~a quality, length ~d" metric length))))))

(defun quality-string (quality metric)
  "Wrapper for ENCODE-QUALITY that encodes QUALITY, an array of bytes
representing base quality scores, as a string using an encoder
appropriate for METRIC, a quality metric ( :SANGER , :SOLEXA
or :ILLUMINA )."
  (let ((encoder (ecase metric
                   (:sanger #'encode-phred-quality)
                   (:solexa #'encode-solexa-quality)
                   (:illumina #'encode-illumina-quality))))
    (encode-quality quality encoder)))

(defun encode-quality (quality encoder)
  "Encodes QUALITY, an array of bytes representing base quality
scores, as a string using function ENCODER."
  (declare (optimize (speed 3)))
  (declare (type (simple-array quality-score (*)) quality)
           (type function encoder))
  (map-into (make-string (length quality) :element-type 'base-char)
            encoder quality))

(defun decode-quality (quality decoder)
  "Decodes the QUALITY, a string, as into a new array using function
DECODER."
  (declare (optimize (speed 3)))
  (declare (type simple-string quality)
           (type function decoder))
  (map-into (make-array (length quality) :element-type 'quality-score)
            decoder quality))

(defun ensure-encoded (vector encoder element-type)
  "If VECTOR is a string, returns an encoded token vector of element
type (unsigned-byte 4) of the same length, otherwise returns
VECTOR. ENCODER is the encoding function used to convert characters to
\(unsigned-byte 4\)."
  (if (stringp vector)
      (when (plusp (length vector))
        (map-into (make-array (length vector) :element-type element-type)
                  encoder vector))
    vector))

(defun ensure-decoded-quality (quality metric)
  "If QUALITY is a string, returns a decoded quality vector of
integers of the same length, otherwise returns QUALITY. METRIC is
either :SANGER, :SOLEXA or :ILLUMINA , denoting the quality metric to
be used when decoding."
  (if (stringp quality)
      (let ((decoder (ecase metric
                       (:sanger #'decode-phred-quality)
                       (:solexa #'decode-solexa-quality)
                       (:illumina #'decode-illumina-quality))))
        (decode-quality quality decoder))
    quality))

(defun token-subsequence (tokens start end)
  "Returns a subsequence of TOKENS between indices START and END."
  (let* ((length (length tokens))
         (end (or end length)))
    (%check-token-range (length tokens) start end)
    (replace (make-array (- end start)
                         :element-type (array-element-type tokens))
             tokens :start2 start :end2 end)))

(defun complement-tokens (tokens comp-fn &optional (start 0) end)
  "Returns a complemented copy of TOKENS populated with elements
from TOKENS that have been transformed by COMP-FN, starting at the
first element, or index START, and continuing to the last residue, or
index END."
  (let ((end (or end (length tokens))))
    (map-into (make-array (- end start)
                          :element-type (array-element-type tokens))
              comp-fn tokens)))

(defun %to-string-virtual (seq char start end token-case)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (%check-token-range length start end)
      (make-string (- end start) :element-type 'base-char
                   :initial-element (ecase token-case
                                      ((nil) char)
                                      (:lower (char-downcase char))
                                      (:upper char))))))

(defun %to-string-encoded (seq decoder start end token-case)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (let ((str (make-string (- end start) :element-type 'base-char)))
        (when (< 0 (length str))
          (copy-array (slot-value seq 'vector) start (1- end)
                      str 0 decoder))
        (ecase token-case
          ((nil) str)
          (:lower str)
          (:upper (nstring-upcase str)))))))

(declaim (inline %check-token-range))
(defun %check-token-range (length start &optional end)
  (declare (optimize (speed 3)))
  (let ((end (or end length)))
    (declare (type vector-index length start end))
    (check-arguments (<= 0 start end length) (start end)
                     "must satisfy (<= 0 start end ~d)" length)))

(declaim (inline %aref-dna-4bit))
(defun %aref-dna-4bit (vector i)
  (declare (optimize (speed 3)))
  (declare (type (encoded-residues 4) vector))
  (decode-dna-4bit (aref vector i)))

(defun (setf %aref-dna-4bit) (value vector i)
  (declare (optimize (speed 3)))
  (declare (type (encoded-residues 4) vector))
  (setf (aref vector i) (encode-dna-4bit value)))

(declaim (inline %aref-rna-4bit))
(defun %aref-rna-4bit (vector i)
  (declare (optimize (speed 3)))
  (declare (type (encoded-residues 4) vector))
  (decode-rna-4bit (aref vector i)))

(defun (setf %aref-rna-4bit) (value vector i)
  (declare (optimize (speed 3)))
  (declare (type (encoded-residues 4) vector))
  (setf (aref vector i) (encode-rna-4bit value)))

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
  (with-slots (vector)
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

(defun maybe-subsequence (seq start &optional end)
  (with-accessors ((length length-of))
      seq
    (let ((end (or end length)))
      (if (and (zerop start) (= length end))
          seq
        (subsequence seq start end)))))
