
(in-package :bio-sequence)

(defvar *seq-len-print-limit* 50)

(deftype encoded-bio-symbols (n)
  `(simple-array (unsigned-byte ,n) *))

(deftype bio-symbol-subscript ()
  '(and fixnum (integer 0 *)))

(defun make-simple-seq (class str encoder)
  "Returns a new bio-sequence of CLASS composed of the sequence
symbols in the string STR encoded in 2 bits per symbol by the function
ENCODER."
  (declare (optimize (speed 3) (safety 1)))
  (let ((symbol-seq (make-array (length str)
                                :element-type '(unsigned-byte 2))))
    (declare (type (encoded-bio-symbols 2) symbol-seq)
             (type simple-string str)
             (type function encoder))
    (copy-array str 0 (1- (length str)) symbol-seq 0 encoder)
    (make-instance class :symbol-seq symbol-seq)))

(defun make-iupac-seq (class str encoder)
  "Returns a new bio-sequence of CLASS composed of the sequence
symbols in the string STR encoded in 2 bits per symbol by the function
ENCODER."
  (declare (optimize (speed 3) (safety 1)))
  (let ((symbol-seq (make-array (length str)
                                :element-type '(unsigned-byte 4))))
    (declare (type (encoded-bio-symbols 4) symbol-seq)
             (type simple-string str)
             (type function encoder))
    (copy-array str 0 (1- (length str)) symbol-seq 0 encoder)
    (make-instance class :symbol-seq symbol-seq)))

(defun make-dna-seq (str &key ambiguity)
   "Returns a new DNA-SEQUENCE object with residues specified by
string STR. Base ambiguity may be defined with the :AMBIGUITY
key. Accepted values for :AMBGUITY are NIL (no ambiguity, the default)
and :IUPAC (IUPAC ambiguity)."
   (ccase ambiguity
     ((nil) (make-simple-seq 'simple-dna-sequence str #'encode-dna-2bit))
     (:iupac (make-iupac-seq 'iupac-dna-sequence str #'encode-dna-4bit))))

(defun make-rna-seq (str &key ambiguity)
   "Returns a new RNA-SEQUENCE object with residues specified by
string STR. Base ambiguity may be defined with the :AMBIGUITY
key. Accepted values for :AMBGUITY are NIL (no ambiguity, the default)
and :IUPAC (IUPAC ambiguity)."
   (ccase ambiguity
     ((nil) (make-simple-seq 'simple-rna-sequence str #'encode-rna-2bit))
     (:iupac (make-iupac-seq 'iupac-rna-sequence str #'encode-rna-4bit))))

(defmethod length-of ((seq bio-sequence))
  (let ((symbol-seq (symbol-seq-of seq)))
    (length symbol-seq)))

(defmethod residue-of ((seq simple-dna-sequence) (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-symbol seq subscript (encoded-bio-symbols 2)))

(defmethod (setf residue-of) (value (seq simple-dna-sequence)
                              (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-symbol value seq subscript (encoded-bio-symbols 2)))

(defmethod residue-of ((seq simple-rna-sequence) (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-symbol seq subscript (encoded-bio-symbols 2)))

(defmethod (setf residue-of) (value (seq simple-rna-sequence)
                              (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-symbol value seq subscript (encoded-bio-symbols 2)))

(defmethod residue-of ((seq iupac-dna-sequence) (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-symbol seq subscript (encoded-bio-symbols 4)))

(defmethod (setf residue-of) (value (seq iupac-dna-sequence)
                              (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-symbol value seq subscript (encoded-bio-symbols 4)))

(defmethod residue-of ((seq iupac-rna-sequence) (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (decode-symbol seq subscript (encoded-bio-symbols 4)))

(defmethod (setf residue-of) (value (seq iupac-rna-sequence)
                              (subscript fixnum))
  (declare (optimize (speed 3) (safety 1)))
  (encode-symbol value seq subscript (encoded-bio-symbols 4)))

(defmethod to-string ((seq simple-dna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type bio-symbol-subscript start end))
  (decode-symbol-array seq start end (encoded-bio-symbols 2)))

(defmethod to-string ((seq simple-rna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type bio-symbol-subscript start end))
  (decode-symbol-array seq start end (encoded-bio-symbols 2)))

(defmethod to-string ((seq iupac-dna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type bio-symbol-subscript start end))
  (decode-symbol-array seq start end (encoded-bio-symbols 4)))

(defmethod to-string ((seq iupac-rna-sequence) &optional
                      (start 0)
                      (end (length-of seq)))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type bio-symbol-subscript start end))
  (decode-symbol-array seq start end (encoded-bio-symbols 4)))

(defmethod copy-sequence ((seq bio-sequence))
  (let ((symbol-seq (symbol-seq-of seq)))
    (make-instance (class-of seq) :symbol-seq
                   (make-array (length symbol-seq)
                               :element-type (array-element-type symbol-seq)
                               :initial-contents symbol-seq))))

(defmethod print-object ((obj simple-dna-sequence) stream)
  (print-object-aux "SIMPLE-DNA-SEQUENCE" obj stream))

(defmethod print-object ((obj simple-rna-sequence) stream)
  (print-object-aux "SIMPLE-RNA-SEQUENCE" obj stream))

(defmethod print-object ((obj iupac-dna-sequence) stream)
  (print-object-aux "IUPAC-DNA-SEQUENCE" obj stream))

(defmethod print-object ((obj iupac-rna-sequence) stream)
  (print-object-aux "IUPAC-RNA-SEQUENCE" obj stream))


(defmacro decode-symbol (seq subscript symbol-seq-type)
  "Returns the decoded residue symbol from SUBSCRIPT in SEQ. The
residue symbols are encoded as an array of type SYMBOL-SEQ-TYPE in
SEQ."
  (let ((symbol-seq (gensym))
        (decoder (gensym)))
    `(let ((,symbol-seq (symbol-seq-of ,seq))
           (,decoder (decoder-of ,seq)))
       (declare (type ,symbol-seq-type ,symbol-seq)
                (type function ,decoder))
       (funcall ,decoder (aref ,symbol-seq ,subscript)))))

(defmacro encode-symbol (value seq subscript symbol-seq-type)
  "Sets the residue symbol VALUE at SUBSCRIPT in SEQ. The residue
symbols are encoded as an array of type SYMBOL-SEQ-TYPE in SEQ."
  (let ((symbol-seq (gensym))
        (encoder (gensym)))
    `(let ((,symbol-seq (symbol-seq-of ,seq))
           (,encoder (encoder-of ,seq)))
       (declare (type ,symbol-seq-type ,symbol-seq)
                (type function ,encoder))
       (setf (aref ,symbol-seq ,subscript)
             (funcall ,encoder ,value)))))

(defmacro decode-symbol-array (seq start end symbol-seq-type)
  (let ((symbol-seq (gensym))
        (dest (gensym))
        (source-end (gensym))
        (dest-start (gensym))
        (decoder (gensym)))
    `(let ((,symbol-seq (symbol-seq-of ,seq))
           (,dest (make-array (- ,end ,start)
                              :element-type 'base-char))
           (,source-end (1- ,end))
           (,dest-start 0)
           (,decoder (decoder-of seq)))
       (declare (type ,symbol-seq-type ,symbol-seq)
                (type simple-base-string ,dest)
                (type bio-symbol-subscript ,source-end ,dest-start)
                (type function ,decoder))
       (copy-array ,symbol-seq ,start ,source-end
                   ,dest ,dest-start ,decoder)
       ,dest)))

(defun print-object-aux (name obj stream)
  "Helper function for printing bio-sequence objects."
  (let ((len (length-of obj)))
    (if (<= len *seq-len-print-limit*)
        (format stream "<~a ~a>" name (to-string obj))
      (format stream "<~a, length ~d>" name len))))

(defun benchmark (n)
  (declare (optimize (speed 3) (safety 1)))
  (let* ((seq (make-dna-seq n))
         (symbol-seq (symbol-seq-of seq))
         (encoder (encoder-of seq)))
    (declare (type function encoder)
             (type (encoded-bio-symbols 2) symbol-seq))
    (loop for i of-type fixnum from 0 below n
          do (setf (aref symbol-seq i) (funcall encoder #\g)))
    seq))


(defun benchmark2 (n)
  (declare (optimize (speed 3) (safety 1)))
  (let ((seq (make-dna-seq n)))
    (loop for i of-type fixnum from 0 below n
       do (setf (residue-of seq i) #\g))))


;; (defun read-fasta (stream)
;;   "Requests DNA sequence from Fasta stream."
;;   (let ((header (read-line stream nil nil))
;;         (string (make-array 1000 :adjustable t :fill-pointer 0
;;                             :element-type 'base-char)))
;;     (when header
;;       (do ((peek (peek-char nil stream nil 'eof)
;;                  (peek-char nil stream nil 'eof)))
;;           ((or (eq 'eof peek) (char= peek 'a)))
;;         (let ((char (read-char stream)))
;;           (unless (char= 'x char)
;;             (vector-push-extend char data 100))))
;;       (let* ((dna-seq (make-dna-sequence (length data)))
;;              (symbol-seq (symbol-seq-of dna-seq)))
;;         (dotimes (i (length data))
;;           (setf (address-cell symbol-seq i)
;;                 (encode-dna-2bit (aref data i))))
;;         dna-seq))))

(defun fasta-in-stream-p (stream)
  (char= #\> (peek-char nil stream nil 'eof)))