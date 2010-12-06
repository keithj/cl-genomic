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

(defmethod make-seq-input ((stream line-input-stream)
                           (format (eql :fastq))
                           &key (alphabet :dna) (metric :sanger) parser)
  (let ((parser (or parser
                    (make-instance 'quality-sequence-parser :metric metric))))
    (defgenerator
        (more (has-sequence-p stream format))
        (next (read-fastq-sequence stream alphabet parser)))))

(defmethod make-seq-output ((stream stream) (format (eql :fastq))
                            &key token-case metric)
  (lambda (obj)
    (write-fastq-sequence obj stream :token-case token-case :metric metric)))

(defmethod split-sequence-file (filespec (format (eql :fastq))
                                pathname-gen &key (chunk-size 1))
  (with-seq-input (seqi (pathname filespec) :fastq
                        :parser (make-instance 'raw-sequence-parser))
    (split-from-generator seqi #'write-fastq-sequence chunk-size pathname-gen)))

(defmethod has-sequence-p ((stream line-input-stream)
                           (format (eql :fastq)) &key alphabet)
  (declare (ignore alphabet))
  (let ((seq-header (find-line stream #'content-string-p)))
    (cond ((eql :eof seq-header)
           nil)
          (t
           (unwind-protect
                (check-record (fastq-header-p seq-header) nil
                              "the stream contains non-Fastq data ~s"
                              seq-header)
             (push-line stream seq-header))))))

(defmethod read-fastq-sequence ((stream line-input-stream)
                                (alphabet symbol)
                                (parser bio-sequence-parser))
  (flet ((parse-residues (id s)
           (let ((residues (stream-read-line s))
                 (quality-header (stream-read-line s))
                 (quality (stream-read-line s)))
             (check-record (and (stringp residues)
                                (stringp quality-header)
                                (stringp quality)
                                (fastq-quality-header-p quality-header))
                           (pairlis '(:identity :residues :quality)
                                    (list id residues quality))
                           "invalid Fastq record ~s" id)
             residues)))
    (restart-case
        (let ((seq-header (find-line stream #'content-string-p)))
          (cond ((eql :eof seq-header)
                 (values nil nil))
                (t
                 (check-field (fastq-header-p seq-header) nil seq-header
                              "~s is not recognised as as Fastq header"
                              seq-header)
                 (let* ((identity (string-left-trim '(#\@) seq-header))
                        (residues (parse-residues identity stream)))
                   (begin-object parser)
                   (object-alphabet parser alphabet)
                   (object-identity parser identity)
                   (object-residues parser residues)
                   (values (end-object parser) t)))))
      (skip-bio-sequence ()
        :report "Skip this sequence."
        ;; Restart skips on to the next header
        (let ((line (find-line stream #'fastq-header-p)))
          (push-line stream line))
        (values nil t)))))

(defmethod read-fastq-sequence ((stream line-input-stream)
                                (alphabet symbol)
                                (parser quality-parser-mixin))
  (restart-case
      (let ((seq-header (find-line stream #'content-string-p))
            (seq-len 0))
        (cond ((eql :eof seq-header)
               (values nil nil))
              (t
               (check-field (fastq-header-p seq-header) nil seq-header
                            "~s is not recognised as a Fastq header"
                            seq-header)
               (begin-object parser)
               (object-alphabet parser alphabet)
               (object-identity parser (string-left-trim '(#\@) seq-header))
               (loop
                  for line = (stream-read-line stream)
                  while (not (eql :eof line))
                  until (fastq-quality-header-p line) ; discard this line
                  do (progn
                       (object-residues parser line)
                       (incf seq-len (length line))))
               (loop
                  for line = (stream-read-line stream)
                  while (not (eql :eof line))
                  until (and (fastq-header-p line) (= seq-len qual-len))
                  sum (length line) into qual-len
                  do (object-quality parser line)
                  finally (unless (eql :eof line) ; push back the new header
                            (push-line stream line)))
               (values (end-object parser) t))))
    (skip-bio-sequence ()
      :report "Skip this sequence."
      ;; Restart skips on to the next header
      (let ((line (find-line stream #'fastq-header-p)))
        (unless (eql :eof line)
          (push-line stream line)))
      (values nil t))))

(defmethod write-fastq-sequence ((seq dna-quality-sequence) (stream stream)
                                 &key token-case metric)
  (let ((*print-pretty* nil)
        (metric (or metric (metric-of seq))))
    (write-char #\@ stream)
    (write-line (if (anonymousp seq)
                    ""
                    (identity-of seq)) stream)
    (write-line (nadjust-case (coerce-sequence seq 'string) token-case) stream)
    (write-line "+" stream)
    (write-line (quality-string (quality-of seq) metric) stream)))

(defmethod write-fastq-sequence ((alist list) (stream stream)
                                 &key token-case metric)
  (declare (ignore metric))
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

(defmethod write-fastq-sequence (obj filespec &key token-case metric)
  (with-open-file (stream filespec :direction :output)
    (write-fastq-sequence obj stream :token-case token-case :metric metric)))

(declaim (inline fastq-header-p))
(defun fastq-header-p (str)
  "Returns T if STR is a Fastq header (starts with the character '@'),
or NIL otherwise. This function is sensitive to the well-known design
fault in the Fastq format, being that in general it is not possible to
distinguish a record header from a line of quality data where the
first character is '@'."
  (starts-with-char-p str #\@))

(declaim (inline fastq-quality-header-p))
(defun fastq-quality-header-p (str)
 "Returns T if STR is a Fastq header (starts with the character '@'),
or NIL otherwise."
  (starts-with-char-p str #\+))
