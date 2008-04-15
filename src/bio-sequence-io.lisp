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

(defmethod read-bio-sequence :before (stream format &key alphabet
                                             virtualp)
  (unless alphabet
    (error 'invalid-argument-error
           :params 'alphabet
           :args nil
           :text "an alphabet is required"))
  (unless (member virtualp '(t nil))
    (error 'invalid-argument-error
           :params 'virtualp
           :args virtualp
           :text "expected a value of T or NIL")))


(defun make-chunk-pname (file-pname chunk-number)
  "Returns a new pathname for a file chunk based on a file pathname
FILE-PNAME and an integer CHUNK-NUMBER."
  (make-pathname
   :directory (pathname-directory file-pname)
   :name (concatenate 'string (pathname-name file-pname) "."
                      (princ-to-string chunk-number))
   :type (pathname-type file-pname)))


;; (defmethod begin-object (handler))
;; (defmethod object-relation (handler relation value))
;; (defmethod atomic-property (handler property value))
;; (defmethod sequence-property (handler property index value))
;; (defmethod end-object (handler))

;; (defclass handler ()
;;   ())

;; (defclass bio-sequence-handler (handler identity-mixin)
;;   (token-seq :iniform nil
;;              :initarg :token-seq-chunks
;;              :accessor token-seq-chunks-of))

;; (defmethod begin-object ((handler bio-sequence-handler)))
;; (defmethod end-object ((handler bio-sequence-handler)))

;; (defmethod atomic-property ((handler bio-sequence-handler) 
;;                             (property (eql :identity))
;;                             (value string))
;;   (setf (identity-of handler) value))

;; (defmethod atomic-property ((handler bio-sequence-handler) object
;;                             (property (eql :description))
;;                             (value string))
;;   (setf (description-of handler) value))

;; (defmethod sequence-property ((handler bio-sequence-handler)
;;                               (property (eql :token-seq))
;;                               (index fixnum)
;;                               (value vector))
;;   (push (cons index value) (token-seq-chunks-of handler)))

;; (defmethod sequence-property ((handler bio-sequence-handler)
;;                               (property (eql :quality))
;;                               (index fixnum)
;;                               (value vector))
;;   (push (cons index value) (quality-chunks-of handler)))
