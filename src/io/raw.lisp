;;;
;;; Copyright (c) 2009-2011 Keith James. All rights reserved.
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

(declaim (type fixnum *raw-line-width*))
(defparameter *raw-line-width* 50
  "Line width for printing raw sequence files.")

;; A raw sequence file is defined as one that contains only residue
;; tokens and whitespace. Any whitespace is ignored, meaning that a
;; raw file can contain only one sequence.

(defmethod make-seq-input ((stream character-line-input-stream)
                           (format (eql :raw))
                           &key (alphabet :dna) parser virtual)
  (let ((parser (or parser
                    (cond (virtual
                           (make-instance 'virtual-sequence-parser))
                          (t
                           (make-instance 'simple-sequence-parser))))))
    (defgenerator
        (more (has-sequence-p stream format))
        (next (read-raw-sequence stream alphabet parser)))))

(defmethod make-seq-output ((stream stream) (format (eql :raw))
                            &key token-case)
  (lambda (obj)
    (write-raw-sequence obj stream :token-case token-case)))

(defmethod has-sequence-p ((stream character-line-input-stream)
                           (format (eql :raw)) &key alphabet)
  (declare (ignore alphabet))
  (let ((line (find-line stream #'content-string-p)))
    (cond ((eql :eof line)
           nil)
          (t
           (push-line stream line)
           t))))

(defmethod read-raw-sequence ((stream character-line-input-stream)
                              (alphabet symbol)
                              (parser bio-sequence-parser))
  (let ((first-line (find-line stream #'content-string-p)))
    (cond ((eql :eof first-line)
           (values nil nil))
          (t
           (begin-object parser)
           (object-alphabet parser alphabet)
           (loop
              for line = first-line then (stream-read-line stream)
              until (eql :eof line)
              do (let ((tmp (nsubstitute #\Space #\Tab line)))
                   (if (find #\Space tmp :test #'char=)
                       (dolist (x (string-split tmp #\Space
                                                :remove-empty-substrings t))
                         (object-residues parser x))
                     (object-residues parser line)))
              finally (return (values (end-object parser) nil)))))))

(defmethod write-raw-sequence ((seq bio-sequence) stream &key token-case)
  (let ((*print-pretty* nil)
        (len (length-of seq)))
    (loop
       for i from 0 below len by *raw-line-width*
       do (write-line
           (nadjust-case
            (coerce-sequence seq 'string
                             :start i :end (min len (+ i *raw-line-width*)))
            token-case) stream))))

(defmethod write-raw-sequence ((alist list) stream &key token-case)
  (let* ((*print-pretty* nil)
         (residues (let ((str (or (assocdr :residues alist) "")))
                     (nadjust-case str token-case)))
         (len (length residues)))
    (loop
       for i from 0 below len by *raw-line-width*
       do (write-line residues stream
                      :start i
                      :end (min len (+ i *raw-line-width*))))))

(defmethod write-raw-sequence (obj filespec &key token-case)
  (with-open-file (stream filespec :direction :output
                          :if-exists :supersede)
    (write-raw-sequence obj stream :token-case token-case)))
