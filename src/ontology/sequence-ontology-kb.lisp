;;;
;;; Copyright (c) 2010-2011 Keith James. All rights reserved.
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

(defparameter *bio-sequence-instances* (make-hash-table :test 'equal)
  "Maps bio-sequence identity strings to their corresponding
bio-sequence object and PowerLoom logic object.")
(defparameter *logic-instances* (make-hash-table)
  "Maps PowerLoom logic objects representing bio-sequence instances to
their corresponding bio-sequence object.")

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

(defun term-same (term1 term2)
  (equal (get-name term1) (get-name term2)))

(defun term-parents (child)
  "Return a list of the \"parent\" terms of term CHILD. This is
parenthood in the GFF3 sense; meaning the terms which CHILD is a
part_of."
  (retrieve `(all ?p (and (concept ?p)
                          (concept ,child)
                          (part_of ,child ?p))) :realise :nconc))

(defun term-parent-p (parent child)
  "Returns T if term PARENT is a \"parent\" of term CHILD. This is
parenthood in the GFF3 sense; meaning the CHILD is a part_of the
PARENT."
  (ask `(part_of ,child ,parent)))

;; In asserting that a sequence object is-a thing:
;;
;; - We need to assert the object in the context of a PL module
;; - We can't change the identity-of the object
;;
;; Therefore we map the identity-of the object to the asserted instance
;;
;; (identity-of object) => (list object logic-object)
;; logic-object => object

(defgeneric assert-instance (bio-sequence term)
  (:documentation "Asserts that BIO-SEQUENCE is an instance of TERM."))

(defgeneric retract-instance (bio-sequence term)
  (:documentation "Retracts the assrtion that BIO-SEQUENCE is an
instance of TERM."))

(defmethod assert-instance ((seq bs:bio-sequence) (concept string))
  (evaluate `(assert (,(stella-symbol concept)
                       ,(stella-symbol (bs:identity-of seq))))))

(defmethod assert-instance ((seq bs:bio-sequence) concept)
  (let ((identity (stella-symbol (bs:identity-of seq))))
    (evaluate `(assert (,concept ,identity)))))

(defmethod assert-instance :after ((seq bs:bio-sequence) concept)
  (with-accessors ((identity bs:identity-of))
      seq
    (let* ((name (stella-symbol identity))
           (instance (retrieve `(1 ?i (= ,name ?i)) :realise :nconc)))
      (when instance
        (setf (gethash identity *bio-sequence-instances*) (cons seq instance)
              (gethash (car instance) *logic-instances*) seq)))))

(defmethod retract-instance ((identity string) concept)
  (evaluate `(retract (,concept ,(stella-symbol identity)))))

(defmethod retract-instance ((seq bs:bio-sequence) concept)
  (evaluate `(retract (,concept ,(stella-symbol (bs:identity-of seq))))))

(defmethod retract-instance :after ((identity string) concept)
  (remhash identity *bio-sequence-instances*)
  (remhash (logic-instance identity) *logic-instances*))

(defmethod retract-instance :after ((seq bs:bio-sequence) concept)
  (with-accessors ((identity bs:identity-of))
      seq
    (remhash identity *bio-sequence-instances*)
    (remhash (logic-instance identity) *logic-instances*)))

(defun find-instance (identity)
  (first (gethash identity *bio-sequence-instances*)))

(defun logic-instance (identity)
  (second (gethash identity *bio-sequence-instances*))) 

(defmethod instancep ((seq bs:bio-sequence))
  (find-instance (bs:identity-of seq)))

(defmethod instance-module ((seq bs:bio-sequence))
  (let ((instance (find-instance (bs:identity-of seq))))
    (when instance
      (get-home-module (second instance)))))

(defmethod instance-terms ((seq bs:bio-sequence))
  (when (instancep seq)
    (with-accessors ((identity bs:identity-of))
        seq
      (get-types (first (retrieve `(1 (= ,(stella-symbol identity) ?x))
                                  :realise :nconc)) :realise :nconc))))

;; (assert-instance (bs:make-dna "atggcatcgcatgc" :identity "foo") "SO:0000673")
;; (assert-instance (bs:make-dna "atggc" :identity "bar") "SO:0000147")

;; (let ((logic-transcript (logic-instance "foo"))
;;       (logic-exon (logic-instance "bar")))
;;   (evaluate `(assert (part_of ,logic-exon ,logic-transcript))))

;; At REPL:
;; (load-ontology "/home/keith/dev/lisp/cl-genomic.git/ontology/so_2_4_3.plm")
;; (load-ontology "/home/keith/dev/lisp/cl-genomic.git/ontology/so_addenda.plm")
;;
;; Use syntax:
;; (in-syntax *powerloom-readtable*)
;;
;; These are all equivalent
;; (with-module (:sequence-ontology)
;;    (get-concept "SO:0000001"))
;; (with-module (:sequence-ontology)
;;   @SO:0000001)
;;
;; Switch to module
;; (in-module :sequence-ontology)
;;
;; (retrieve '(all ?x (part_of |SO:0000179| ?x)))
;;
;; With the dollar reader syntax enabled:
;; (retrieve '(all ?x (part_of |SO:0000179| ?x)))
;; (retrieve '(all ?x (part_of $SO:0000179 ?x)))
;; (retrieve '(all ?x (part_of $SO:0000179 ?x)) :realise :nconc)
;; (retrieve '(all (?x ?d) (and (part_of $SO:0000179 ?x)
;;                                 (= (documentation ?x) ?d))))
;;
;; Splice in Lisp values with backquote:
;;
;; (let ((name "TSS"))
;;   (retrieve `(1 (and (concept ?c)
;;                      (= (name ?c) ,name))) :realise :nconc))
;;
;;
;; Get tree of sub-terms. Requires concept argument
;; (subrelation-tree @SO:0000001)
;;
;; Apply term-name function to all nodes to get new tree
;; (traverse * #'term-name)
;;
;; Get all part_of SO:0000179
;; (retrieve '(all ?x (part_of $SO:0000179 ?x)))
;;
;; or just
;; (retrieve '(all (part_of $SO:0000179 ?x)))
;;
;; or flattened
;; (retrieve '(all (part_of $SO:0000179 ?x)) :realise :nconc)
;;
;; or as generator function
;; (retrieve '(all (part_of $SO:0000179 ?x)) :realise nil)
