;;;
;;; Copyright (C) 2009-2009 Keith James. All rights reserved.
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

(defmethod end-section ((parser obo-powerloom-parser))
  (with-accessors ((state state-of) (header header-of)
                   (tag-values tag-values-of) (terms terms-of)
                   (typedefs typedefs-of) (instances instances-of))
      parser
    (ecase state
      (header (setf header tag-values))
      (term (merge-tag-values tag-values terms))
      (typedef (merge-tag-values tag-values typedefs))
      (instance (merge-tag-values tag-values instances))
      ((nil) nil)) ; case for multiple or trailing empty lines
    (setf state nil
          tag-values ())))

(defun convert-obo-powerloom (in-filespec out-filespec module &key
                              (base-concept 'thing) (parent-module "PL-USER"))
  (with-open-file (out out-filespec :direction :output
                       :if-exists :supersede)
    (with-li-stream (in in-filespec)
      (let* ((*default-base-concept* base-concept)
             (parser (read-obo-stream in (make-instance
                                           'obo-powerloom-parser))))
        (write-powerloom-module module parent-module :stream out)
        (write-powerloom-deffunction base-concept 'name 'string out)
        (write-powerloom-relations (typedefs-of parser) :stream out)
        (write-powerloom-concepts (terms-of parser) :stream out)))))

(defun write-powerloom-concept (id &optional (stream *standard-output*))
  (format stream "(defconcept ~a)~%" id))

(defun write-powerloom-doc (id doc &optional (stream *standard-output*))
  (format stream "(assert (documentation ~a ~s))~%" id doc))

(defun write-powerloom-deffunction (id name type
                                    &optional (stream *standard-output*))
  (format stream "(deffunction ~a ((?x ~a)) :->~%   (?~a ~a))~%"
          name id name type))

(defun write-powerloom-subset (id super
                               &optional (stream *standard-output*))
  (format stream "(assert (subset-of ~a ~a))~%" id super))

(defun write-powerloom-assert (form &optional (stream *standard-output*))
  (format stream "(assert ~a)~%" form))

(defun write-powerloom-relation (relation subject object
                                 &optional (stream *standard-output*))
  (format stream "(defrelation ~a ((?x ~a) (?y ~a)))~%"
          relation subject object))

(defun write-powerloom-subrelation (sub-relation super-relation
                                    subject object
                                    &optional (stream *standard-output*))
  (format stream "(defrelation ~a ((?x ~a) (?y ~a))~%   :=> (~a ?x ?y))~%"
          sub-relation subject object super-relation))

(defun write-powerloom-module (module parent-module
                               &key (stream *standard-output*))
  (format stream "(defmodule ~s~%  :includes (~s))~%"
          (string-upcase module) (string-upcase parent-module))
  (format stream "(in-module ~s)~%" (string-upcase module))
  (format stream "(clear-module ~s)~%~%" (string-upcase module)))

(defun write-powerloom-concepts (terms &key (stream *standard-output*))
  (let ((deferred ()))
    (loop
       for id being the hash-keys of terms
       using (hash-value tag-values)
       do (progn
            (write-powerloom-concept id stream)
            (loop
               for (tag . value) in tag-values
               when (string= "intersection_of" tag) collect
               value into intersect
               do (cond ((string= "def" tag)
                         (write-powerloom-doc id value stream))
                        ((string= "name" tag)
                         (write-powerloom-assert
                          (format nil "(= (name ~a) ~s)" id value) stream))
                        ((string= "is_a" tag)
                         (write-powerloom-subset id value stream))
                        ((string= "relationship" tag)
                         (write-powerloom-assert
                          (format nil "(~a ~a ~a)"
                                  (first value) id (second value)) stream))
                        (t
                         nil))
               ;; FIXME - raise parse error when single intersection clause
               finally (when intersect
                         (push (with-output-to-string (s)
                                 (write-powerloom-intersect id intersect s)
                                 (terpri s))
                               deferred)))
            (terpri stream)))
    (dolist (d deferred)
      (write-line d stream))
    (terpri stream)))

(defun write-powerloom-intersect (id intersect
                                  &optional (stream *standard-output*))
  (format stream (txt "(assert (forall ?term~%~10t(=>"
                      "(and~%~{~a~^~%~})~%~12t(~a ?term))))")
          (mapcar (lambda (x)
                    (if (atom x)
                        (format nil "~16t(subset-of ?term ~a)" x)
                      (format nil "~16t(~a ?term ~a)" (first x) (second x))))
                  intersect) id))

(defun write-powerloom-relations (typedefs &key
                                  (base-concept *default-base-concept*)
                                  (stream *standard-output*))
  (loop
     for id being the hash-keys of typedefs
     using (hash-value tag-values)
     do (progn
          (unless (assocdr "is_a" tag-values :test #'string=)
            (write-powerloom-relation id base-concept base-concept stream))
          (loop
             for (tag . value) in tag-values
             do (cond ((string= "def" tag)
                       (write-powerloom-doc id value stream))
                       ((string= "is_a" tag)
                        (write-powerloom-subrelation id value
                                                     base-concept base-concept
                                                     stream))
                      ((string= "is_transitive" tag)
                       (when (string= "true" value)
                         (write-powerloom-assert
                          (format nil "(transitive ~a)" id) stream)))
                      ((string= "is_symmetric" tag)
                       (when (string= "true" value)
                         (write-powerloom-assert
                          (format nil "(symmetric ~a)" id) stream)))
                      (t
                       nil)))
          (terpri stream))))

(defun merge-tag-values (tag-values table)
  "Removes the id tag-value cons from the list TAG-VALUES and merges
the remainder into a list, together with any previous data for that
id. The list is inserted into TABLE using id as the key."
  (let ((id (assocdr "id" tag-values :test #'string=))
        (others (loop
                   for (tag . value) in (mapcar #'read-value tag-values)
                   unless (string= "id" tag)
                   collect (cons tag value))))
    (multiple-value-bind (val presentp)
        (gethash id table)
      (if presentp
          (setf (gethash id table) (nconc others val))
        (setf (gethash id table) others)))
    table))
