;;;
;;; Copyright (C) 2008-2009 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
;;
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

(in-package :stella)

(defmodule "GENE-ONTOLOGY"
  :includes ("PL-USER"))

(in-module "GENE-ONTOLOGY")
(clear-module "GENE-ONTOLOGY")

(in-dialect KIF)

(defconcept molecular-function ((?t obo-term))
  :<=> (and (obo-term ?t)
            (= (namespace ?t) "molecular_function")))

(deffunction name ((?t obo-term)) :-> (?s STRING)
             :documentation "The term name.")

(deffunction namespace ((?t obo-term)) :-> (?n STRING)
             :documentation "The term namespace.")

;; (retrieve all (molecular-function ?x))
;; (retrieve all (subset-of |GO:0008150| ?x))

;; GO:0003674 molecular_function

(in-package :bio-sequence)

(defun convert-ontology (obo-filespec powerloom-filespec)
  (with-open-file (out powerloom-filespec :direction :output
                       :element-type 'character
                       :if-exists :supersede)
    (prin1 '(in-module "GENE-ONTOLOGY") out)
    (terpri out)
    (prin1 '(defconcept obo-term) out)
    (terpri out)
    (prin1 '(deffunction name ((?t obo-term)) :-> (?n STRING)
             :documentation "The term name.") out)
    (terpri out)
    (prin1 '(deffunction namespace ((?t obo-term)) :-> (?n STRING)
             :documentation "The term namespace.") out)
    (terpri out)
    (prin1 '(assert (subset-of |GO:0008150| obo-term)) out)
    (terpri out)
    (with-ascii-li-stream (stream obo-filespec)
      (let ((id nil))
        (loop
           as line = (stream-read-line stream)
           while (not (eql :eof line))
           do (cond ((id-tag-p line)
                     (setf id (intern (parse-tag line)))
                     (prin1 `(defconcept ,id) out)
                     (terpri out))
                    ((name-tag-p line)
                     (prin1 `(assert (= (name ,id)
                                      ,(subseq line 6))) out)
                     (terpri out))
                    ((namespace-tag-p line)
                     (prin1 `(assert (= (namespace ,id)
                                      ,(subseq line 11))) out)
                     (terpri out))
;;                     ((def-tag-p line)
;;                      (prin1 `(assert (documentation ,id
;;                                       ,(subseq line 5))) out)
;;                      (terpri out))
                    ((isa-tag-p line)
                     (let ((isa (intern (parse-tag line))))
                       (prin1 `(assert (subset-of ,isa ,id)) out)
                       (terpri out)))
                    (t
                     nil)))))))

(defun term-stanza-p (str)
  (starts-with-string-p str "[Term]"))

(defun id-tag-p (str)
  (starts-with-string-p str "id:"))

(defun name-tag-p (str)
  (starts-with-string-p str "name:"))

(defun namespace-tag-p (str)
  (starts-with-string-p str "namespace:"))

(defun isa-tag-p (str)
  (starts-with-string-p str "is_a:"))

(defun def-tag-p (str)
  (starts-with-string-p str "def:"))

(defun relationship-tag-p (str)
  (starts-with-string-p str "relationship:"))

(defun parse-tag (str)
  (second (string-split #\Space str)))


(defparameter *default-module* pli:null)
(defparameter *default-environment* pli:null)

(defmacro with-environment ((environment) &body body)
  `(let ((*default-environment* ,environment))
    ,@body))

(defmacro with-module ((module) &body body)
  `(let ((*default-module* (find-module ,module)))
    ,@body))

(defun find-module (symbol &optional (environment *default-environment*))
  (pli:get-module (symbol-name symbol) environment))

(defun convert-symbol (symbol &key (module *default-module*)
                       (environment *default-environment*))
  (pli:create-symbol (symbol-name symbol) module environment))

(defun convert-tree (tree &key (module *default-module*)
                     (environment *default-environment*))
  (cond ((symbolp tree)
         (convert-symbol tree :module module :environment environment))
        ((listp tree)
         (mapcar (lambda (x)
                   (convert-tree x :module module :environment environment))
                 tree))
        (t
         tree)))

(defun retrieve (query &key (module *default-module*)
                 (environment *default-environment*))
  (pli:retrieve (convert-tree query :module module :environment environment)
                module environment))

(retrieve '(all (subset-of |GO:0008150| ?x))
          :module (find-module :gene-ontology))

(let ((x *))
  (retrieve `((= (name ,x) ?x)) :module (find-module :gene-ontology)))

(retrieve '(all (?y ?n) (and (subset-of ?y |SO:0000336|) (= (name ?y) ?n)))
          :module (find-module :sequence-ontology))
