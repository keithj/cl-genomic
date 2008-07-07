;;;
;;; Copyright (C) 2007-2008, Keith James. All rights reserved.
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

(define-condition bio-sequence-io-error (io-error bio-sequence-error)
  ((text :initform nil
         :initarg :text
         :reader text-of
         :documentation "Error message text."))
  (:report (lambda (condition stream)
             (format stream "IO error~@[: ~a~]"
                     (text-of condition))))
  (:documentation "An error that is raised when performing stream of
file IO on bio-sequences."))

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

