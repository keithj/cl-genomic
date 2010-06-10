;;;
;;; Copyright (C) 2010 Keith James. All rights reserved.
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

(in-package :bio-ontology)

(in-syntax *powerloom-readtable*)

(defun traverse (tree fn) ; move this to utilities
  (cond ((null tree)
         nil)
        ((atom tree)
         (funcall fn tree))
        (t
         (cons (traverse (first tree) fn)
               (traverse (rest tree) fn)))))

;;; term-* functions take a concept argument (an OBO term).
;;;
;;; find-* functions take a string argument that is used to find a
;;; concept (an OBO term).

(defun term-name (concept)
  "Returns the name of term CONCEPT."
  (first (retrieve `(1 (and (concept ,concept)
                            (= (name ,concept) ?n))) :realise :nconc)))

(defun term-doc (concept)
  "Returns the documentation of term CONCEPT."
  (first (retrieve `(1 (and (concept ,concept)
                            (= (documentation ,concept) ?d))) :realise :nconc)))

(defun find-term (name)
  "Returns the term named NAME, or NIL if there os not such term."
  (first (retrieve `(1 (and (concept ?c)
                            (= (name ?c) ,name))) :realise :nconc)))

(defun find-doc (name)
  "Returns the documentation of a term with NAME, or NIL if there os
not such term."
  (first (retrieve `(1 (?d ?c)
                       (and (concept ?c)
                            (= (name ?c) ,name)
                            (= (documentation ?c) ?d))) :realise :nconc)))

(defun term-parents (child)
  "Return a list of the \"parent\" terms of term CHILD. This is
parenthood in the GFF3 sense; meaning the terms which CHILD is a
part_of."
  (retrieve `(all ?p (and (concept ?p)
                          (concept ,child)
                          (part_of ,child ?p))) :realise :nconc))

(defun term-parent-p (parent child)
  (ask `(part_of ,child ,parent)))

;; At REPL:
;; (load-ontology "/home/keith/dev/lisp/cl-genomic.git/ontology/sofa_2_4_2.plm")
;; (load-ontology "/home/keith/dev/lisp/cl-genomic.git/ontology/sofa_addenda.plm")
;;
;; Use syntax:
;; (in-syntax *powerloom-readtable*)
;;
;; These are all equivalent
;; (with-module (:sofa)
;;    (get-concept "SO:0000001"))
;; (with-module (:sofa)
;;   $"SO:0000001")
;; (with-module (:sofa)
;;   $|SO:0000001|)
;;
;; Switch to module
;; (in-module :sofa)
;;
;; With the dollar reader syntax enabled:
;; (retrieve '(all ?x (part_of $"SO:0000179" ?x)))
;; (retrieve '(all ?x (part_of $|SO:0000179| ?x)))
;; (retrieve '(all ?x (part_of $|SO:0000179| ?x)) :realise :nconc)
;; (retrieve '(all (?x ?n ?d) (and (part_of $"SO:0000179" ?x)
;;                                 (= (documentation ?x) ?d))))
;;
;; Splice in Lisp values with backquote:
;;
;; (let ((name "TSS"))
;;   (retrieve `(1 (and (concept ?c)
;;                      (= (name ?c) ,name))) :realise :nconc))
;;
;;
;; Get tree of sub-terms
;; (subrelation-tree $|SO:0000001|)
;;
;; Apply term-name function to all nodes to get new tree
;; (traverse * #'term-name)
;;
;; Get all part_of SO:0000179
;; (retrieve '(all ?x (part_of $|SO:0000179| ?x)))
;;
;; or just
;; (retrieve '(all (part_of $|SO:0000179| ?x)))
;;
;; or flattened
;; (retrieve '(all (part_of $|SO:0000179| ?x)) :realise :nconc)
;;
;; or as generator function
;; (retrieve '(all (part_of $|SO:0000179| ?x)) :realise nil)
