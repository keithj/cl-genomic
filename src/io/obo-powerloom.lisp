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

;;; Parser method for collecting OBO data
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

;;; This conversion has not been well tested. The output was
;;; syntactically correct (loaded into PowerLoom without errors) for
;;; the Sequence Ontology release 2.4. The semantics of the result
;;; could probably be improved.
(defun convert-obo-powerloom (obo-filespec plm-filespec module &key
                              (base-concept 'thing) (parent-module :pl-user)
                              (if-exists :supersede))
  "Converts OBO ontology data to PowerLoom declarations.

Arguments:

- obo-filespec (pathname designator): The OBO file to be read.
- plm-filespec (pathname designator): The PowerLoom file to be written.
- module (symbol): A symbol naming the new PowerLoom module into which
  the declarations will be written.

Key:

- base-concept (symbol): A symbol naming the base PowerLoom concept
  that will be used as the range of all converted OBO
  relations. Defaults to 'thing .
- parent-module (symbol): A symbol naming the parent PowerLoom module
  of MODULE. Defaults to :pl-user .

- if-exists (symbol): One of NIL , :error or :supersede indicating the
  behaviour if PLM-FILESPEC already exists.

Returns:

- A filespec  (pathname designator) for the PowerLoom file."
  (with-open-file (plm plm-filespec :direction :output
                       :if-exists if-exists)
    (with-li-stream (obo obo-filespec)
      (let* ((*default-base-concept* base-concept)
             (parser (read-obo-stream obo (make-instance
                                           'obo-powerloom-parser)))
             (tmpl (txt ";;; Created ~d-~2,'0d-~2,'0d ~d:~2,'0d:~2,'0d"
                        "from OBO file ~a by"
                        "cl-genomic OBO to PowerLoom converter.~%")))
        (multiple-value-bind
              (sec min hour date month year)
            (get-decoded-time)
          (write-line (format nil tmpl year month date
                              hour min sec (parse-file obo-filespec)) plm))
        (write-powerloom-module module parent-module :stream plm)
        (write-powerloom-deffunction base-concept 'name 'string plm)
        (write-powerloom-relations (typedefs-of parser) :stream plm)
        (write-powerloom-concepts (terms-of parser) :stream plm))))
  plm-filespec)

(defun write-powerloom-concept (id &optional (stream *standard-output*))
  "Writes a defconcept declaration for OBO term ID to STREAM."
  (format stream "(defconcept ~a)~%" id))

(defun write-powerloom-doc (id doc &optional (stream *standard-output*))
  "Writes a documentation declaration DOC for OBO term ID to STREAM."
  (format stream "(assert (documentation ~a ~s))~%" id doc))

(defun write-powerloom-deffunction (id name type
                                    &optional (stream *standard-output*))
  "Writes a deffunction declaration NAME for OBO term ID to STREAM."
  (format stream "(deffunction ~a ((?x ~a)) :->~%   (?~a ~a))~%"
          name id name type))

(defun write-powerloom-subset (id super
                               &optional (stream *standard-output*))
  "Writes a subset-of declaration asserting OBO term ID to be a subset
of OBO term SUPER to STREAM."
  (format stream "(assert (subset-of ~a ~a))~%" id super))

(defun write-powerloom-assert (form &optional (stream *standard-output*))
  "Writes a PowerLoom declaration asserting FORM to STREAM."
  (format stream "(assert ~a)~%" form))

(defun write-powerloom-relation (relation subject object
                                 &optional (stream *standard-output*))
  "Writes a PowerLoom defrelation declaration declaring RELATION
between OBO terms SUBJECT and OBJECT to STREAM."
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
  "Writes PowerLoom concepts using OBO term data from hash-table TERMS
to STREAM. The hash-table keys are OBO term IDs, while the hash-values
are conses containing OBO tags and values."
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
