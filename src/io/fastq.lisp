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

(defmethod read-bio-sequence-alist ((stream binary-line-input-stream)
                                    (format (eql :fastq))
                                    &key alphabet ambiguity virtualp
                                    (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line stream #'byte-fastq-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (seq quality-header quality)
            (read-fastq-record stream #'byte-fastq-quality-header-p)
          (declare (ignore quality-header))
          (let ((alist (make-quality-alist
                        (make-sb-string seq-header 1)
                        alphabet ambiguity
                        :token-seq (unless virtualp (make-sb-string seq))
                        :length (when virtualp (length seq))
                        :quality (make-sb-string quality))))
            (if callback-supplied-p
                (apply callback alist callback-args)
              alist)))
      nil)))

(defmethod read-bio-sequence-alist ((stream character-line-input-stream)
                                    (format (eql :fastq))
                                    &key alphabet ambiguity virtualp
                                    (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line stream #'char-fastq-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (seq quality-header quality)
            (read-fastq-record stream #'char-fastq-quality-header-p)
          (declare (ignore quality-header))
          (let ((alist (make-quality-alist
                        (string-left-trim '(#\@) seq-header)
                        alphabet ambiguity
                        :token-seq (unless virtualp seq)
                        :length (when virtualp (length seq))
                        :quality quality)))
            (if callback-supplied-p
                (apply callback alist callback-args)
              alist)))
      nil)))

(defmethod read-bio-sequence (stream (format (eql :fastq))
                              &key alphabet ambiguity virtualp metric)
  (read-bio-sequence-alist stream format :alphabet alphabet
                           :ambiguity ambiguity :virtualp virtualp
                           :callback #'make-quality-seq-fastq
                           :callback-args (list metric)))


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

(defun make-quality-seq-fastq (alist metric)
  "Callback which accepts a ALIST and creates a new
dna-quality-sequence with quality METRIC."
  (make-quality-seq :alphabet (assocdr :alphabet alist)
                    :ambiguity (assocdr :ambiguity alist)
                    :token-seq (assocdr :token-seq alist)
                    :quality (assocdr :quality alist)
                    :identity (assocdr :identity alist)
                    :metric metric))

(defun write-alist-fastq (alist &optional output-stream)
  "Callback which accepts an ALIST and writes it to OUTPUT-STREAM as a
Fastq format record. OUTPUT-STREAM defaults to *standard-output*."
  (let ((*print-pretty* nil))
    (write-char #\@ output-stream)
    (write-line (assocdr :identity alist) output-stream)
    (write-line (assocdr :token-seq alist) output-stream)
    (write-line "+" output-stream)
    (write-line (assocdr :quality alist) output-stream))
  t)

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

;; (defun write-n-fastq (stream n chunk-pname)
;;   "Reads up to N Fastq records from STREAM and writes them into a new
;; file of pathname CHUNK-PNAME. Returns the number of records actually
;; written, which may be 0 if the stream contained to further records."
;;   (if (more-lines-p stream)
;;       (with-open-file (out chunk-pname :direction :output
;;                        :if-exists :supersede
;;                        :element-type 'base-char)
;;         (do ((fq (read-bio-sequence-alist stream :fastq :alphabet :dna
;;                                           :callback #'write-alist-fastq
;;                                           :callback-args (list out))
;;                  (read-bio-sequence-alist stream :fastq :alphabet :dna
;;                                           :callback #'write-alist-fastq
;;                                           :callback-args (list out)))
;;              (count 1 (1+ count)))
;;             ((or (null fq)
;;                  (= count n)) count)))
;;     0))

;; (defun write-n-fastq (stream n chunk-pname)
;;   "Reads up to N Fastq records from STREAM and writes them into a new
;; file of pathname CHUNK-PNAME. Returns the number of records actually
;; written, which may be 0 if the stream contained to further records."
;;   (with-open-file (out chunk-pname :direction :output
;;                    :if-exists :supersede
;;                    :element-type 'base-char)
;;     (loop
;;        for count from 0 below n
;;        with written = 0
;;        do (let ((fq (read-bio-sequence-alist stream :fastq :alphabet :dna)))
;;             (if (null fq)
;;                 (return written)
;;               (progn
;;                 (write-alist-fastq fq out)
;;                 (incf written))))
;;        finally (return written))))

(defun write-n-fastq (stream n chunk-pname)
  "Reads up to N Fastq records from STREAM and writes them into a new
file of pathname CHUNK-PNAME. Returns the number of records actually
written, which may be 0 if the stream contained to further records."
  (let ((written 
         (with-open-file (out chunk-pname :direction :output
                          :if-exists :supersede
                          :element-type 'base-char)
           (loop
              for count from 0 below n
              for fq = (read-bio-sequence-alist stream :fastq :alphabet :dna)
              then (read-bio-sequence-alist stream :fastq :alphabet :dna)
              while fq
              do (write-alist-fastq fq out)
              finally (return count)))))
    (when (zerop written)
      (delete-file chunk-pname))
    written))

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
