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

(defmethod read-bio-sequence ((stream line-input-stream)
                              (format (eql :fastq))
                              &key alphabet virtualp metric)
  (read-seq-datum stream format :alphabet alphabet
                  :virtualp virtualp
                  :callback #'make-quality-seq-from-datum
                  :callback-args (list metric)))

(defmethod write-bio-sequence ((seq dna-quality-sequence)
                               stream (format (eql :fastq))  &key token-case)
  (let ((*print-pretty* nil)
        (encoder (ecase (metric-of seq)
                   (:phred #'encode-phred-quality)
                   (:illumina #'encode-illumina-quality))))
    (write-char #\@ stream)
    (write-line (identity-of seq) stream)
    (write-line (to-string seq :token-case token-case) stream)
    (write-line "+" stream)
    (write-line (encode-quality (quality-of seq) encoder) stream)))

(defmethod read-seq-datum ((stream binary-line-input-stream)
                           (format (eql :fastq))
                           &key alphabet virtualp
                           (callback nil callback-supplied-p)
                           callback-args)
  (let ((seq-header (find-line stream #'byte-fastq-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (seq quality-header quality)
            (read-fastq-record stream #'byte-fastq-quality-header-p)
          (declare (ignore quality-header))
          (let ((datum (make-quality-datum
                        (make-sb-string seq-header 1)
                        alphabet
                        :token-seq (unless virtualp (make-sb-string seq))
                        :length (when virtualp (length seq))
                        :quality (make-sb-string quality))))
            (if callback-supplied-p
                (apply callback datum callback-args)
              datum)))
      nil)))

(defmethod read-seq-datum ((stream character-line-input-stream)
                           (format (eql :fastq))
                           &key alphabet virtualp
                           (callback nil callback-supplied-p)
                           callback-args)
  (let ((seq-header (find-line stream #'char-fastq-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (seq quality-header quality)
            (read-fastq-record stream #'char-fastq-quality-header-p)
          (declare (ignore quality-header))
          (let ((datum (make-quality-datum
                        (string-left-trim '(#\@) seq-header)
                        alphabet
                        :token-seq (unless virtualp seq)
                        :length (when virtualp (length seq))
                        :quality quality)))
            (if callback-supplied-p
                (apply callback datum callback-args)
              datum)))
      nil)))

(defmethod write-seq-datum ((stream stream) (format (eql :fastq)) datum)
  (let ((*print-pretty* nil))
    (write-char #\@ stream)
    (write-line (seq-datum-identity datum) stream)
    (write-line (seq-datum-token-seq datum) stream)
    (write-line "+" stream)
    (write-line (seq-datum-quality datum) stream))
  t)

(defmethod filter-seq-datum ((stream line-input-stream)
                             (format (eql :fastq)) pred out)
  "Reads Fastq records from STREAM and writes only those for which
PRED returns T to stream character stream OUT. PRED should be a
function that accepts a single argument of a standard seq-datum and
returns T if the read should be removed, or NIL otherwise."
  (do ((fq (read-seq-datum stream :fastq :alphabet :dna)
           (read-seq-datum stream :fastq :alphabet :dna))
       (num-written 0))
      ((null fq) num-written)
    (when (not (funcall pred fq))
      (write-seq-datum out :fastq fq)
      (incf num-written))))

(defun read-fastq-record (stream qual-header-validate-fn)
  "Reads the body of a Fastq record that follows the header from
STREAM, validates the quality header with the predicate
QUAL-HEADER-VALIDATE-FN and returns three vector values: the sequence,
the quality header and the quality."
  (let ((seq (stream-read-line stream))
        (quality-header (stream-read-line stream))
        (quality (stream-read-line stream)))
    (unless (and (vectorp seq)
                 (funcall qual-header-validate-fn quality-header)
                 (vectorp quality))
      (error 'malformed-record-error :text
             "Incomplete Fastq record."))
    (values seq quality-header quality)))


(defun split-fastq-file (filespec chunk-size)
  "Splits Fastq file identified by FILESPEC into automatically named
chunks, each, except the last file, containing up to CHUNK-SIZE
records."
  (let ((file-pname (pathname filespec)))
    (with-open-file (in file-pname :direction :input
                     :element-type '(unsigned-byte 8))
      (do* ((stream (make-line-input-stream in))
            (chunk-count 0 (1+ chunk-count))
            (chunk-pname (make-chunk-pname file-pname chunk-count)
                         (make-chunk-pname file-pname chunk-count))
            (n (write-n-fastq stream chunk-size chunk-pname)
               (write-n-fastq stream chunk-size chunk-pname)))
           ((zerop n))))))

(defun write-n-fastq (stream n chunk-pname)
  "Reads up to N Fastq records from STREAM and writes them into a new
file of pathname CHUNK-PNAME. Returns the number of records actually
written, which may be 0 if the stream contained to further records."
  (let ((num-written 
         (with-open-file (out chunk-pname :direction :output
                          :if-exists :supersede
                          :element-type 'base-char)
           (loop
              for count from 0 below n
              for fq = (read-seq-datum stream :fastq :alphabet :dna)
              then (read-seq-datum stream :fastq :alphabet :dna)
              while fq
              do (write-seq-datum out :fastq fq)
              finally (return count)))))
    (when (zerop num-written)
      (delete-file chunk-pname))
    num-written))

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
