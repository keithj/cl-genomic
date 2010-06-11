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

(in-package :cl-user)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun powerloom-symbols ()
    "Returns a list of the symbols required for the PowerLoom Lisp
API."
    '(#:clear-module
      #:get-current-module
      #:get-home-module
      #:get-child-modules
      #:get-parent-modules

      #:get-name
      #:get-domain
      #:get-range
      #:get-arity

      #:is-default
      #:is-enumerated-collection
      #:is-enumerated-list
      #:is-enumerated-set
      #:is-false
      #:is-float
      #:is-integer
      #:is-logic-object
      #:is-number
      #:is-string
      #:is-true
      #:is-unknown)))

(defmacro define-bio-ontology-package ()
  `(defpackage :bio-ontology
     (:use #:common-lisp #:deoxybyte-utilities)
     (:nicknames #:bo)
     (:import-from #:pli ,@(powerloom-symbols))
     (:export

      #:load-ontology
      #:in-syntax

      #:with-environment

      #:with-module
      #:in-module
      #:modulep
      #:get-module

      #:get-concept

      #:is-subrelation
      #:is-a
      #:is-true-binary-proposition
      #:is-true-proposition
      #:is-true-unary-proposition

      #:evaluate
      #:ask
      #:retrieve

      #:subrelation-tree
      #:traverse
   
      #:term-name
      #:term-doc
      #:term-parents
      #:term-parent-p

      #:find-term
      #:find-doc)))

(define-bio-ontology-package)
