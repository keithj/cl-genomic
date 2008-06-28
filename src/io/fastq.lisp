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

(defmethod make-input-gen ((stream line-input-stream)
                           (format (eql :fastq))
                           &key (alphabet :dna) (metric :phred) parser)
  (let* ((parser (or parser
                    (make-instance 'quality-sequence-parser
                                   :metric metric)))
         (current (read-fastq-sequence stream alphabet parser)))
    (lambda (op)
      (ecase op
        (:current current)
        (:next (prog1
                   current
                 (setf current
                       (read-fastq-sequence stream alphabet parser))))
        (:more (not (null current)))))))

(defmethod make-output-con ((stream stream) (format (eql :fastq))
                            &key token-case)
  (lambda (bio-sequence)
    (write-fastq-sequence bio-sequence stream :token-case token-case)))


(defmethod split-sequence-file (filespec (format (eql :fastq))
                                generator &key (chunk-size 1))
  (let ((file-pathname (pathname filespec)))
    (with-open-file (stream file-pathname :direction :input
                     :element-type 'base-char
                     :external-format :ascii)
      (split-sequence-stream stream #'write-n-fastq chunk-size generator))))

(defmethod read-fastq-sequence ((stream binary-line-input-stream)
                                (alphabet symbol)
                                (parser quality-parser-mixin))
  (let ((seq-header (find-line stream #'byte-fastq-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (residues quality-header quality)
            (parse-fastq-record stream #'byte-fastq-quality-header-p)
          (declare (ignore quality-header))
          (begin-object parser)
          (object-alphabet parser alphabet)
          (object-identity parser (make-sb-string seq-header 1))
          (object-residues parser (make-sb-string residues))
          (object-quality parser (make-sb-string quality))
          (end-object parser))
      nil)))

(defmethod read-fastq-sequence ((stream character-line-input-stream)
                                (alphabet symbol)
                                (parser quality-parser-mixin))
  (let ((seq-header (find-line stream #'char-fastq-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (residues quality-header quality)
            (parse-fastq-record stream #'char-fastq-quality-header-p)
          (declare (ignore quality-header))
          (begin-object parser)
          (object-alphabet parser alphabet)
          (object-identity parser (string-left-trim '(#\@) seq-header))
          (object-residues parser residues)
          (object-quality parser quality)
          (end-object parser))
      nil)))

(defmethod write-fastq-sequence ((seq dna-quality-sequence) stream
                                 &key token-case)
  (let ((*print-pretty* nil)
        (encoder (ecase (metric-of seq)
                   (:phred #'encode-phred-quality)
                   (:illumina #'encode-illumina-quality))))
    (write-char #\@ stream)
    (write-line (identity-of seq) stream)
    (write-line (to-string seq :token-case token-case) stream)
    (write-line "+" stream)
    (write-line (encode-quality (quality-of seq) encoder) stream)))

(defun write-raw-fastq (raw stream)
  "Writes sequence data RAW to STREAM in Fastq format. The alist RAW
must contain keys and values as created by {defclass raw-sequence-parser} ."
  (let ((*print-pretty* nil))
    (write-char #\@ stream)
    (write-line (assocdr :identity raw) stream)
    (write-line (assocdr :residues raw) stream)
    (write-line "+" stream)
    (write-line (assocdr :quality raw) stream)))

(defun parse-fastq-record (stream qual-header-validate-fn)
  "Reads the body of a Fastq record that follows the header from
STREAM, validates the quality header with the predicate
QUAL-HEADER-VALIDATE-FN and returns three vector values: the sequence,
the quality header and the quality."
  (let ((residues (stream-read-line stream))
        (quality-header (stream-read-line stream))
        (quality (stream-read-line stream)))
    (unless (and (vectorp residues)
                 (funcall qual-header-validate-fn quality-header)
                 (vectorp quality))
      (error 'malformed-record-error :text
             "Incomplete Fastq record."))
    (values residues quality-header quality)))

(defun split-fastq-file (filespec chunk-size generator)
  (split-sequence-file filespec :fastq generator :chunk-size chunk-size))

(defun write-n-fastq (stream n pathname)
  "Reads up to N Fastq records from STREAM and writes them into a new
file of PATHNAME. Returns the number of records actually written,
which may be 0 if STREAM contained no further records."
  (let ((gen (make-input-gen stream :fastq :alphabet :dna
                             :parser (make-instance
                                      'raw-sequence-parser))))
    (write-n-raw-sequences gen #'write-raw-fastq n pathname)))

(defun byte-fastq-header-p (bytes)
  "Returns T if BYTES are a Fastq header (start with the character
code for '@'), or NIL otherwise."
  (starts-with-byte-p bytes (char-code #\@)))

(defun char-fastq-header-p (str)
  "Returns T if STR is a Fastq header (starts with the character '@'),
or NIL otherwise."
  (starts-with-char-p str #\@))

(defun byte-fastq-quality-header-p (bytes)
  "Returns T if BYTES are a Fastq quality header (start with the
character code for '+'), or NIL otherwise."
  (starts-with-byte-p bytes (char-code #\+)))

(defun char-fastq-quality-header-p (str)
 "Returns T if STR is a Fastq header (starts with the character '@'),
or NIL otherwise."
  (starts-with-char-p str #\+))

(defun concat-quality-arrays (quality-arrays)
  (let ((new-quality (make-array (reduce #'+ quality-arrays :key #'length)
                                 :element-type 'quality-score))
        (num-arrays (length quality-arrays)))
    (do ((i 0 (1+ i))
         (offset 0))
        ((= i num-arrays) new-quality)
      (let ((quality-array (aref quality-arrays i)))
        (unless (zerop (length quality-array))
          (copy-array quality-array 0 (1- (length quality-array))
                      new-quality offset)
          (incf offset (length quality-array)))))))
