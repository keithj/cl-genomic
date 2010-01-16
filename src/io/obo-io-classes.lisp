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

(defclass obo-parser ()
  ((state :initform 'header
          :accessor state-of
          :documentation "The parser state tracks the type of document
section from which the tags and values were derived.")
   (tag-values :initform ()
               :accessor tag-values-of
               :documentation "The tags and values for a header or
single stanza."))
  (:documentation "A parser for processing Open Biomedical Ontologies
format version 1.2"))

(defclass obo-powerloom-parser (obo-parser)
  ((header :initform ()
           :accessor header-of
           :documentation "The header tag-values.")
   (terms :initform (make-hash-table :test #'equal)
          :accessor terms-of
          :documentation "All OBO terms, indexed by id.")
   (typedefs :initform (make-hash-table :test #'equal)
             :accessor typedefs-of
             :documentation "All OBO typedefs, indexed by id.")
   (instances :initform (make-hash-table :test #'equal)
              :accessor instances-of
              :documentation "All OBO instances, indexed by id."))
  (:documentation "An OBO parser that collects state for conversion
  into PowerLoom format."))
