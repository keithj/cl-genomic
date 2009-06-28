;;;
;;; Copyright (C) 2007-2009 Keith James. All rights reserved.
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

(define-condition bio-sequence-error (error)
  ()
  (:documentation "The parent type of all bio-sequence error conditions."))

(define-condition bio-sequence-warning (warning)
  ()
  (:documentation "The parent type of all bio-sequence warning
conditions."))

;; (define-condition bio-sequence-io-error (io-error
;;                                          bio-sequence-error)
;;   ((text :initform nil
;;          :initarg :text
;;          :reader text-of
;;          :documentation "Error message text."))
;;   (:report (lambda (condition stream)
;;              (format stream "IO error~@[: ~a~]"
;;                      (text-of condition))))
;;   (:documentation "An error that is raised when performing stream IO
;; on bio-sequences."))

;; (define-condition bio-sequence-parse-error (general-parse-error
;;                                             bio-sequence-io-error)
;;   ((text :initform nil
;;          :initarg :text
;;          :reader text-of
;;          :documentation "Error message text."))
;;   (:report (lambda (condition stream)
;;              (format stream "IO error~@[: ~a~]"
;;                      (text-of condition))))
;;   (:documentation "An error that is raised when performing stream IO
;; on bio-sequences."))

(define-condition bio-sequence-op-error (invalid-operation-error
                                         bio-sequence-error)
  ((text :initform nil
         :initarg :text
         :reader text-of
         :documentation "Error message text."))
  (:report (lambda (condition stream)
             (format stream "Invalid bio-sequence operation~@[: ~a~]"
                     (text-of condition))))
  (:documentation "An error that is raised when performing operations
on bio-sequences."))

(define-condition initiator-codon-error (bio-sequence-op-error)
  ((codon :initform nil
          :initarg :codon
          :reader codon-of
          :documentation "The invalid codon.")
   (genetic-code :initform nil
                 :initarg :genetic-code
                 :reader genetic-code-of
                 :documentation "The genetic code used to translate."))
  (:report (lambda (condition stream)
             (format stream "Codon ~a is not an initiator in ~a"
                     (codon-of condition)
                     (genetic-code-of condition))))
  (:documentation "An error that is raised when attempting to translate
a non-initiator codon as an initiator."))

(define-condition translation-error (bio-sequence-op-error)
  ((seq :initform nil
        :initarg :sequence
        :reader sequence-of
        :documentation "The translated sequence.")
   (start :initform nil
          :initarg :start
          :reader start-of
          :documentation "The translation start position.")
   (end :initform nil
        :initarg :end
        :reader end-of
        :documentation "The translation end position.")
   (genetic-code :initform nil
                 :initarg :genetic-code
                 :reader genetic-code-of
                 :documentation "The genetic code used to translate."))
  (:report (lambda (condition stream)
             (format stream "Invalid translation of ~a~@[: ~a~]"
                     (sequence-of condition)
                     (text-of condition))))
  (:documentation "An error that is raised when attempting an invalid
  translation of a sequence."))
