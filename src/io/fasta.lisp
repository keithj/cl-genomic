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

(defparameter *fasta-line-width* 50
  "Line width for printing Fasta files.")

(defparameter *token-cache-extend* 256
  "The number of elements by which the token cache is extended when it
becomes full of chunks of sequence tokens.")


(defmethod read-bio-sequence-datum ((stream binary-line-input-stream)
                                    (format (eql :fasta))
                                    &key alphabet ambiguity virtualp
                                    (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line stream #'byte-fasta-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (identity description)
            (parse-fasta-header (make-sb-string seq-header))
          (let ((datum (make-seq-datum identity alphabet ambiguity
                                       :token-seq (unless virtualp
                                                    (concat-fasta-chunks
                                                     stream
                                                     #'byte-fasta-header-p
                                                     #'concat-into-sb-string))
                                       :length (when virtualp
                                                 (count-fasta-residues
                                                  stream
                                                  #'byte-fasta-header-p))
                                       :description description)))
            (if callback-supplied-p
                (apply callback datum callback-args)
              datum)))
      nil)))

(defmethod read-bio-sequence-datum ((stream character-line-input-stream)
                                    (format (eql :fasta))
                                    &key alphabet ambiguity virtualp
                                    (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line stream #'char-fasta-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (identity description)
            (parse-fasta-header seq-header)
          (let ((datum (make-seq-datum identity alphabet ambiguity
                                      :token-seq (unless virtualp
                                                    (concat-fasta-chunks
                                                     stream
                                                     #'char-fasta-header-p
                                                     #'concat-strings))
                                       :length (when virtualp
                                                 (count-fasta-residues
                                                  stream
                                                  #'char-fasta-header-p))
                                       :description description)))
            (if callback-supplied-p
                (apply callback datum callback-args)
              datum)))
      nil)))

(defmethod read-bio-sequence (stream (format (eql :fasta))
                              &key alphabet ambiguity virtualp)
  (read-bio-sequence-datum stream format :alphabet alphabet
                           :ambiguity ambiguity :virtualp virtualp
                           :callback #'make-seq-from-datum))

(defun count-fasta-residues (stream header-p-fn)
  "Returns an integer which is the number of Fasta sequence residues
read from line-input-stream STREAM. Lines are read until the next
Fasta header is encountered (detected by HEADER-P-FN) or until the end
of the stream is reached."
  (let ((num-residues 0))
    (loop
       as line = (stream-read-line stream)
       while (vectorp line)
       until (funcall header-p-fn line)
       do (incf num-residues (length line))
       finally (when (vectorp line)
                 (push-line stream line))) ; push back the new header
    num-residues))

(defun concat-fasta-chunks (stream header-p-fn concat-fn)
  "Returns a string of concatenated Fasta sequence chunks read from
line-input-stream STREAM. Lines are read until the next Fasta header
is encountered (detected by HEADER-P-FN) or until the end of the
stream is reached. The lines are concatenated using CONCAT-FN."
  (let ((seq-cache (make-array 0 :adjustable t :fill-pointer t)))
    (loop
       as line = (stream-read-line stream)
       and cache-extend = (max 256 (floor (/ (length seq-cache) 2)))
       while (vectorp line)
       until (funcall header-p-fn line)
       do (vector-push-extend line seq-cache cache-extend)
       finally (when (vectorp line)
                 (push-line stream line))) ; push back the new header
    (when (zerop (length seq-cache))
      (error 'malformed-record-error :text "Incomplete Fasta record."))
    (funcall concat-fn seq-cache)))

(defun write-datum-fasta (datum &optional output-stream)
  "A callback which writes sequence data that has been parsed into an
ALIST to OUTPUT-STREAM in Fasta format."
  (let ((description (seq-datum-description datum))
        (*print-pretty* nil))
    (write-char #\> output-stream)
    (if (zerop (length description))
        (write-line (seq-datum-identity datum) output-stream)
      (progn
        (write-string (seq-datum-identity datum) output-stream)
        (write-char #\Space output-stream)
        (write-line description output-stream)))
    (write-wrapped-string (seq-datum-token-seq datum)
                          *fasta-line-width* output-stream))
  t)

(defun byte-fasta-header-p (bytes)
  "Returns T if BYTES are a Fasta header (start with the character
code for '>'), or NIL otherwise."
  (starts-with-byte-p bytes (char-code #\>)))

(defun char-fasta-header-p (str)
  "Returns T if STR is a Fasta header (starts with the character
'>'), or NIL otherwise."
  (starts-with-char-p str #\>))

(defun parse-fasta-header (str)
  "Performs a basic parse of a Fasta header string STR by removing the
leading '>' character and splitting the line on the first space(s)
into identity and description. This function supports pathological
cases where the identity, description, or both are empty strings."
  (multiple-value-bind (split index)
      (split-sequence #\Space str :count 1 :remove-empty-subseqs t)
    (let ((str-len (length str))
          (str-elt-type (array-element-type str))
          (identity-str (car split)))
      (values
       (if (> (length identity-str) 1)
           (adjust-array identity-str (- (length identity-str) 1)
                         :displaced-to identity-str
                         :displaced-index-offset 1)
         (make-string 0 :element-type str-elt-type))
       (if (< index str-len)
           (adjust-array str (- str-len index)
                         :displaced-to str
                         :displaced-index-offset index)
         (make-string 0 :element-type str-elt-type))))))
