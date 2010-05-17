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

(in-package :bio-sequence)

(in-syntax *powerloom-readtable*)

;; (load-file "/home/keith/dev/lisp/cl-genomic.git/ontology/sofa_2_4.plm")

;; These are all equivalent
;; (with-module (:sofa)
;;    (get-concept "SO:0000001"))
;; (with-module (:sofa)
;;   $"SO:0000001")
;; (with-module (:sofa)
;;   $|SO:0000001|)

(defun traverse (tree fn)
  (cond ((null tree)
         nil)
        ((atom tree)
         (funcall fn tree))
        (t
         (cons (traverse (first tree) fn)
               (traverse (rest tree) fn)))))

(defun term-name (concept)
  (first (next-tuple
          (retrieve `(1 ?name (and (concept ,concept)
                                   (= (name ,concept) ?name)))))))

(defun term-doc (concept)
  (first (next-tuple
          (retrieve `(1 ?doc (and (concept ,concept)
                                  (= (documentation ,concept) ?doc)))))))

;; At REPL:
;;
;; Use syntax
;; (in-syntax *powerloom-readtable*)
;; Switch to module
;; (in-module :sofa)
;; Get tree of sub-terms
;; (subrelation-tree $|SO:0000001|)
;; Apply term-name function to all nodes to get new tree
;; (traverse * #'term-name)

(defun find-by-name (name)
  (first (next-tuple
          (retrieve `(1 ?concept (and (concept ?concept)
                                      (= (name ?concept) ,name)))))))

(defun find-all-by-name (name)
  (nconc-tuples (retrieve `(all (?concept ?doc)
                                (and (concept ?concept)
                                     (= (name ?concept) ,name)
                                     (= (documentation ?concept) ?doc))))))

;; Get all part_of SO:0000179
;; (collect-tuples (retrieve `(all ?x  (part_of ,$|SO:0000179| ?x))))
