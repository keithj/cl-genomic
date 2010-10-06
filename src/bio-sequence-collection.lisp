;;;
;;; Copyright (C) 2009-2010 Keith James. All rights reserved.
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

;; Need to group sequences that are related in some way:
;; From same organism
;; From same sequencing experiment (may be multiple organisms)
;; From same bioinformatic experiment
;; Other

;; Use cases:
;;
;; A GFF3 file
;; A genome sequence
;; A lane of short-read sequencing data

;; Just throwing around some ideas here; this code is not for use


;; get sequence by key (id?)

;; get sequence by predicate (similar to Lisp find fn) - in fact, a
;; find-type fn covers finding by key and other attributes.


;; Also implement
;; (size-of coll)
;; (memberp seq coll :test #'test)
;; 

(defgeneric find-sequence (seq collection &key key))

(defgeneric add-sequence (seq collection &key key))

(defgeneric remove-sequence (seq collection &key key))

(defclass bio-sequence-collection ()
  ((sequences :initform nil
              :initarg :sequences
              :reader sequences-of)
   (alphabet :initarg :alphabet
             :reader alphabet-of)))

(defclass bio-sequence-vector (bio-sequence-collection)
  ((sequences :initform (make-array 1024 :adjustable t :fill-pointer 0))))

(defclass bio-sequence-table (bio-sequence-collection)
  ((sequences :initform (make-hash-table :test #'equal))))

(defmethod find-sequence (seq (collection bio-sequence-vector) &key key)
  (with-slots (sequences test)
      collection
    (find seq sequences :key key :test #'equal)))

(defmethod find-sequence (seq (collection bio-sequence-table) &key key)
  (with-slots (sequences)
      collection
    (gethash (funcall-if-fn key seq) collection)))


(defun index-sequence-file (filespec format alphabet)
  (let ((index (merge-pathnames (make-pathname :type "index") filespec))
        (data (merge-pathnames (make-pathname :type "data") filespec)))
    (with-open-file (data-stream data :direction :output
                                 :if-exists :supersede
                                 :element-type 'base-char
                                 :external-format :ascii)
      (with-open-file (index-stream index :direction :output
                                 :if-exists :supersede
                                 :element-type 'base-char
                                 :external-format :ascii)
        (let ((parser (make-instance 'indexing-sequence-parser
                                     :stream data-stream)))
          (with-seq-input (seqi input-stream format :alphabet alphabet
                                :parser parser)
            (loop
               for seq = (next seqi)
               while (has-more-p seqi)
               finally (prin1 (index-of parser) index-stream))))))))

(defun index-sequence-file2 (filespec format alphabet)
  (let ((index (make-instance 'tc:tc-hdb))
        (index-file (merge-pathnames (make-pathname :type "index") filespec))
        (data (merge-pathnames (make-pathname :type "data") filespec)))
    (tc:dbm-open index (namestring index-file) :write :create)
    (tc:dbm-optimize index :leaf 512 :non-leaf 256 :bucket-size 100000000
                     :rec-align 4 :free-pool 10 :opts '(:large))
    (tc::dbm-xmsize index 536870912)
    (unwind-protect
         (with-open-file (data-stream data :direction :output
                                      :if-exists :supersede
                                      :element-type 'base-char
                                      :external-format :ascii)
           (let ((parser (make-instance 'indexing-sequence-parser
                                        :stream data-stream)))
             (with-seq-input (seqi input-stream format :alphabet alphabet
                                   :parser parser)
               (loop
                  for seq = (next seqi)
                  for i = 0 then (1+ i)
                  with time = (get-universal-time)
                  with total-time = 0
                  while (has-more-p seqi)
                  do (progn
                       (tc:dbm-put index (identity-of seq) (princ-to-string i))
                       (when (zerop (rem i 100000))
                         (let* ((now (get-universal-time))
                                (interval (- now time)))
                           (setf time now)
                           (incf total-time interval)
                           (format t "~d records in ~d seconds, total ~d seconds~%"
                                   i interval total-time)))))))))
    (tc:dbm-close index)))

(defun index-sequence-file3 (filespec format alphabet)
  (let ((index (merge-pathnames (make-pathname :type "index") filespec))
        (data (merge-pathnames (make-pathname :type "data") filespec)))
    (with-open-file (data-stream data :direction :output
                                 :if-exists :supersede
                                 :element-type 'base-char
                                 :external-format :ascii)
      (with-open-file (index-stream index :direction :output
                                    :if-exists :supersede
                                    :element-type 'base-char
                                    :external-format :ascii)
        (let ((parser (make-instance 'indexing-sequence-parser
                                     :stream data-stream)))
          (with-seq-input (seqi input-stream format :alphabet alphabet
                                :parser parser)
            (loop
               for rlen = (- (offset-of parser) (parsed-length-of parser))
               for seq = (next seqi)
               while (has-more-p seqi)
               do (prog1
                      (princ (identity-of seq) index-stream)
                    (write-char #\Tab index-stream)
                    (princ rlen index-stream)
                    (terpri index-stream)))))))))


(defun lookup-sequence (identity filespec)
  (let ((index (with-open-file (s (merge-pathnames
                                   (make-pathname :type "index") filespec))
                 (read s)))
        (data (merge-pathnames (make-pathname :type "data") filespec)))
    (let ((entry (binary-search index identity :test #'string<
                                               :key #'first)))
      (dxn:with-mapped-vector (vector 'dxn:mapped-vector-char
                                      :filespec data
                                      :length (reduce #'+ index
                                                      :key #'third))
        (loop
           for i from (second entry) below (+ (second entry)
                                              (third entry))
           collect (code-char (dxn:mref vector i)) into bases
           finally (return (values identity bases)))))))
