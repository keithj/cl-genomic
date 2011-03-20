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

(declaim (type array-index *pure-buffer-length*))
(defparameter *pure-buffer-length* 4096)

;; A pure sequence file is defined as one that contains only residue
;; tokens and an optional new-line character at the end of the file. A
;; pure sequence file can contain only a single sequence and is
;; suitable for mmap operations.

(defmethod make-seq-input ((stream character-line-input-stream)
                           (format (eql :pure))
                           &key (alphabet :dna) parser virtual)
  (let ((parser (or parser
                    (cond (virtual
                           (make-instance 'virtual-sequence-parser))
                          (t
                           (make-instance 'simple-sequence-parser))))))
    (defgenerator
        (more (has-sequence-p stream format))
        (next (read-pure-sequence stream alphabet parser)))))

(defmethod make-seq-output ((stream stream) (format (eql :pure))
                            &key token-case)
  (lambda (obj)
    (write-pure-sequence obj stream :token-case token-case)))

(defmethod has-sequence-p ((stream character-line-input-stream)
                           (format (eql :pure)) &key alphabet)
  (declare (ignore alphabet))
  (more-lines-p stream))

(defmethod read-pure-sequence ((stream character-line-input-stream)
                               (alphabet symbol)
                               (parser bio-sequence-parser))
  (cond ((more-lines-p stream)
         (begin-object parser)
         (object-alphabet parser alphabet)
         (loop
            with buffer = (make-array *pure-buffer-length*
                                      :element-type 'base-char)
            for n = (stream-read-sequence stream buffer 0 *pure-buffer-length*)
            while (plusp n)
            do (if (< n *pure-buffer-length*) ; Expect newline at end
                   (object-residues parser (string-right-trim
                                            '(#\Newline) (subseq buffer 0 n)))
                 (object-residues parser buffer))
            finally (return (values (end-object parser) nil))))
        (t
         (values nil nil))))

(defmethod write-pure-sequence ((seq bio-sequence) (stream stream)
                                &key token-case)
  (let ((*print-pretty* nil))
    (write-string (nadjust-case (coerce-sequence seq 'string) token-case)
                  stream)))

(defmethod write-pure-sequence ((alist list) (stream stream) &key token-case)
  (let* ((*print-pretty* nil)
         (residues (let ((str (or (assocdr :residues alist) "")))
                     (nadjust-case str token-case))))
    (write-string residues stream)))

(defmethod write-pure-sequence (obj filespec &key token-case)
  (with-open-file (stream filespec :direction :output
                          :if-exists :supersede)
    (write-pure-sequence obj stream :token-case token-case)))
