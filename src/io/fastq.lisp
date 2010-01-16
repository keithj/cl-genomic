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

(defmethod make-seq-input ((stream character-line-input-stream)
                           (format (eql :fastq))
                           &key (alphabet :dna) (metric :phred) parser)
  (let ((parser (or parser
                    (make-instance 'quality-sequence-parser
                                   :metric metric))))
    (defgenerator
        :next (read-fastq-sequence stream alphabet parser)
        :more (has-sequence-p stream format))))

(defmethod make-seq-output ((stream stream) (format (eql :fastq))
                            &key token-case)
  (lambda (obj)
    (write-fastq-sequence obj stream :token-case token-case)))

(defmethod split-sequence-file (filespec (format (eql :fastq))
                                pathname-gen &key (chunk-size 1))
  (let ((file-pathname (pathname filespec)))
    (with-ascii-li-stream (stream file-pathname)
      (split-from-generator
       (make-seq-input stream :fastq
                       :parser (make-instance 'raw-sequence-parser))
       #'write-fastq-sequence chunk-size pathname-gen))))

(defmethod has-sequence-p ((stream character-line-input-stream)
                           (format (eql :fastq)) &key alphabet)
  (declare (ignore alphabet))
  (let ((seq-header (find-line stream #'content-string-p)))
    (cond ((eql :eof seq-header)
           nil)
          ((char-fastq-header-p seq-header)
           (push-line stream seq-header)
           t)
          (t
           (push-line stream seq-header)
           (error 'malformed-record-error
                  :record seq-header
                  :text "the stream does not contain Fastq data")))))

(defmethod read-fastq-sequence ((stream character-line-input-stream)
                                (alphabet symbol)
                                (parser bio-sequence-parser))
  (restart-case
      (let ((seq-header (find-line stream #'content-string-p)))
        (cond ((eql :eof seq-header)
               (values nil nil))
              ((char-fastq-header-p seq-header)
               (multiple-value-bind (residues quality-header quality)
                   (parse-fastq-record stream #'char-fastq-quality-header-p)
                 (declare (ignore quality-header quality))
                 (begin-object parser)
                 (object-alphabet parser alphabet)
                 (object-identity parser (string-left-trim '(#\@) seq-header))
                 (object-residues parser residues)
                 (values (end-object parser) t)))
              (t
               (error 'malformed-record-error
                      :text (format nil
                                    "~s is not recognised as as Fastq header"
                                    seq-header)))))
    (skip-bio-sequence ()
      :report "Skip this sequence."
      ;; Restart skips on to the next header
      (let ((line (find-line stream #'char-fastq-header-p)))
        (push-line stream line))
      (values nil t))))

(defmethod read-fastq-sequence ((stream character-line-input-stream)
                                (alphabet symbol)
                                (parser quality-parser-mixin))
  (restart-case
      (let ((seq-header (find-line stream #'content-string-p)))
        (cond ((eql :eof seq-header)
               (values nil nil))
              ((char-fastq-header-p seq-header)
               (multiple-value-bind (residues quality-header quality)
                   (parse-fastq-record stream #'char-fastq-quality-header-p)
                 (declare (ignore quality-header))
                 (begin-object parser)
                 (object-alphabet parser alphabet)
                 (object-identity parser (string-left-trim '(#\@) seq-header))
                 (object-residues parser residues)
                 (object-quality parser quality)
                 (values (end-object parser) t)))
              (t
               (error 'malformed-record-error
                      :text (format nil
                                    "~s is not recognised as as Fastq header"
                                    seq-header)))))
    (skip-bio-sequence ()
      :report "Skip this sequence."
      ;; Restart skips on to the next header
      (let ((line (find-line stream #'char-fastq-header-p)))
        (push-line stream line))
      (values nil t))))

(defmethod write-fastq-sequence ((seq dna-quality-sequence) (stream stream)
                                 &key token-case)
  (let ((*print-pretty* nil))
    (write-char #\@ stream)
    (write-line (if (anonymousp seq)
                    ""
                  (identity-of seq)) stream)
    (write-line (nadjust-case (coerce-sequence seq 'string) token-case) stream)
    (write-line "+" stream)
    (write-line (quality-string (quality-of seq) (metric-of seq)) stream)))

(defmethod write-fastq-sequence ((alist list) (stream stream) &key token-case)
  (let ((*print-pretty* nil)
        (residues (let ((str (or (assocdr :residues alist) "")))
                    (nadjust-case str token-case)))
        (quality (or (assocdr :quality alist) ""))
        (identity (or (assocdr :identity alist) "")))
    (write-char #\@ stream)
    (write-line identity stream)
    (write-line residues stream)
    (write-line "+" stream)
    (write-line quality stream)))

(defmethod write-fastq-sequence (obj filespec &key token-case)
  (with-open-file (stream filespec :direction :output)
    (write-fastq-sequence obj stream :token-case token-case)))

(defun write-fastq-alist (alist stream)
  "Writes sequence data ALIST to STREAM in Fastq format. ALIST
must contain keys and values as created by {defclass raw-sequence-parser} ."
  (let ((*print-pretty* nil))11
    (write-char #\@ stream)
    (write-line (or (assocdr :identity alist) "") stream)
    (write-line (or (assocdr :residues alist) "") stream)
    (write-line "+" stream)
    (write-line (or (assocdr :quality alist) "") stream)))

(defun parse-fastq-record (stream qual-header-validate-fn)
  "Reads the body of a Fastq record that follows the header from
STREAM, validates the quality header with the predicate
QUAL-HEADER-VALIDATE-FN and returns three vector values: the sequence,
the quality header and the quality."
  (declare (optimize (speed 3)))
  (declare (type function qual-header-validate-fn))
  (let ((residues (stream-read-line stream))
        (quality-header (stream-read-line stream))
        (quality (stream-read-line stream)))
    (unless (and (vectorp residues)
                 (funcall qual-header-validate-fn quality-header)
                 (vectorp quality))
      (error 'malformed-record-error :text "Incomplete Fastq record."))
    (values residues quality-header quality)))

(defun char-fastq-header-p (str)
  "Returns T if STR is a Fastq header (starts with the character '@'),
or NIL otherwise. This function is sensitive to the well-known design
fault in the Fastq format, being that in general it is not possible to
distinguish a record header from a line of quality data where the
first character is '@'."
  (starts-with-char-p str #\@))

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
