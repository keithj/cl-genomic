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

(defmethod read-bio-sequence-alist ((stream binary-line-input-stream)
                                    (format (eql :fasta))
                                    alphabet ambiguity
                                    &optional (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line stream #'byte-fasta-header-p)))
    (if seq-header
        (multiple-value-bind (identity description)
            (parse-fasta-header (make-sb-string seq-header))
          (let ((alist (make-seq-alist
                        identity alphabet ambiguity
                        (concat-fasta-chunks stream
                                             #'byte-fasta-header-p
                                             #'concat-into-sb-string)
                        description)))
            (if callback-supplied-p
                (apply callback alist callback-args)
              alist)))
      nil)))

(defmethod read-bio-sequence-alist ((stream character-line-input-stream)
                                    (format (eql :fasta))
                                    alphabet ambiguity
                                    &optional (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line stream #'char-fasta-header-p)))
    (if seq-header
        (multiple-value-bind (identity description)
            (parse-fasta-header seq-header)
          (let ((alist (make-seq-alist
                        identity alphabet ambiguity
                        (concat-fasta-chunks stream
                                             #'char-fasta-header-p
                                             #'concat-strings)
                        description)))
            (if callback-supplied-p
                (apply callback alist callback-args)
              alist)))
      nil)))


(defun concat-fasta-chunks (stream header-p-fn concat-fn)
  (let ((seq-cache (make-array 0 :adjustable t :fill-pointer t)))
    (loop
       as line = (stream-read-line stream)
       and cache-extend = (max 256 (floor (/ (length seq-cache) 2)))
       while line
       until (funcall header-p-fn line)
       do (vector-push-extend line seq-cache cache-extend)
       finally (when line
                 (push-line stream line)))
    (when (zerop (length seq-cache))
      (error 'malformed-record-error :text "Incomplete Fasta record."))
    (funcall concat-fn seq-cache)))


(defun write-alist-fasta (alist &optional output-stream)
  "A callback which writes sequence data that has been parsed into an
ALIST to OUTPUT-STREAM in Fasta format."
  (let ((description (assocdr :description alist))
        (*print-pretty* nil))
    (write-char #\> output-stream)
    (if (zerop (length description))
        (write-line (assocdr :identity alist) output-stream)
      (progn
        (write-string (assocdr :identity alist) output-stream)
        (write-char #\Space output-stream)
        (write-line description output-stream)))
    (write-wrapped-string (assocdr :token-seq alist)
                          *fasta-line-width* output-stream))
  t)

(defun byte-fasta-header-p (bytes)
  (starts-with-byte-p bytes (char-code #\>)))

(defun char-fasta-header-p (str)
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
