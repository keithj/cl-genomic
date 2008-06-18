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


(defmethod make-input-gen ((stream line-input-stream)
                           (format (eql :fasta))
                           &key (alphabet :dna) parser virtual)
  (let* ((parser (or parser
                     (cond (virtual
                            (make-instance 'virtual-sequence-parser))
                           (t
                            (make-instance 'simple-sequence-parser)))))
         (current (read-fasta-sequence stream alphabet parser)))
    (lambda (op)
        (ecase op
          (:current current)
          (:next (prog1
                     current
                   (setf current
                         (read-fasta-sequence stream alphabet parser))))
          (:more (not (null current)))))))

(defmethod make-output-con ((stream stream) (format (eql :fasta))
                            &key token-case)
  (lambda (bio-sequence)
    (write-fasta-sequence bio-sequence stream :token-case token-case)))

(defmethod read-fasta-sequence ((stream binary-line-input-stream)
                                (alphabet symbol)
                                (parser bio-sequence-parser))
  (let ((seq-header (find-line stream #'byte-fasta-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (identity description)
            (parse-fasta-header (make-sb-string seq-header))
          (begin-object parser)
          (object-alphabet parser alphabet)
          (object-identity parser identity)
          (object-description parser description)
          (loop
             as line = (stream-read-line stream)
             with offset = 0
             while (vectorp line)
             until (byte-fasta-header-p line)
             do (progn
                  (object-residues parser (make-sb-string line))
                  (incf offset (length line)))
             finally (when (vectorp line) ; push back the new header
                       (push-line stream line)))
          (end-object parser))
      nil)))

(defmethod read-fasta-sequence ((stream character-line-input-stream)
                                (alphabet symbol)
                                (parser bio-sequence-parser))
  (let ((seq-header (find-line stream #'char-fasta-header-p)))
    (if (vectorp seq-header)
        (multiple-value-bind (identity description)
            (parse-fasta-header seq-header)
          (begin-object parser)
          (object-alphabet parser alphabet)
          (object-identity parser identity)
          (object-description parser description)
          (loop
             as line = (stream-read-line stream)
             with offset = 0
             while (vectorp line)
             until (char-fasta-header-p line)
             do (progn
                  (object-residues parser line)
                  (incf offset (length line)))
             finally (when (vectorp line) ; push back the new header
                       (push-line stream line)))
          (end-object parser))
      nil)))

(defmethod write-fasta-sequence ((seq bio-sequence) stream
                                 &key token-case) 
  (let ((*print-pretty* nil)
        (len (length-of seq)))
    (write-char #\> stream)
    (write-line (identity-of seq) stream)
    (loop
       for i from 0 below len by *fasta-line-width*
       do (write-line
           (to-string seq
                      :start i :end (min len (+ i *fasta-line-width*))
                      :token-case token-case)
           stream))))

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
