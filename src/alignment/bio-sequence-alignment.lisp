;;;
;;; Copyright (C) 2008 Keith James. All rights reserved.
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

(defclass aligned-mixin ()
  ((aligned :initarg :aligned
            :reader aligned-of
            :documentation "An aligned copy of a sequence."))
  (:documentation "A mixin describing an aligned copy of a
  sequence."))

(defclass na-alignment-interval (na-sequence-interval aligned-mixin)
  ())

(defclass aa-alignment-interval (aa-sequence-interval aligned-mixin)
  ())

(defclass alignment ()
  ((intervals :initform nil
               :initarg :intervals
               :accessor intervals-of
               :documentation "A list of aligned bio-sequence
               intervals."))
  (:documentation "A sequence alignment."))

(defgeneric aligned-length-of (aligned)
  (:documentation "Returns the aligned length of a sequence or
  sequence interval. This value may include the length contributed by
  gaps."))

;;; Initialization methods
(defmethod initialize-instance :after ((interval na-alignment-interval) &key)
  (with-slots (lower upper reference aligned)
      interval
    (when reference
      (unless (= (- (length-of interval)
                    (num-gaps-of reference :start lower :end upper))
                 (- (length-of aligned)
                    (num-gaps-of aligned)))
        (error 'invalid-argument-error
               :params '(reference lower upper aligned)
               :args (list reference lower upper aligned)
               :text "the reference subsequence and aligned sequence were different lengths, excluding gaps")))))

;;; Printing methods
(defmethod print-object ((interval na-alignment-interval) stream)
  (with-slots (lower upper strand)
      interval
    (format stream "#<NA-ALIGNMENT-INTERVAL ~a ~a ~a" lower upper strand)))

(defmethod print-object ((interval aa-alignment-interval) stream)
  (with-slots (lower upper)
      interval
    (format stream "#<AA-ALIGNMENT-INTERVAL ~a ~a " lower upper)))

(defmethod print-object ((alignment alignment) stream)
  (with-slots (intervals)
      alignment
    (write-line "#<ALIGNMENT" stream)
    (dolist (interval intervals)
      (with-slots (lower upper aligned)
          interval
        (format stream "~7d ~a ~7a~%" lower (to-string aligned) upper)))
    (write-line ">" stream)))

;;; Implementation methods
(defmethod aligned-length-of ((aligned aligned-mixin))
  (with-slots (aligned)
      aligned
    (length-of aligned)))

(defmethod to-string ((interval na-alignment-interval) &key
                      (start 0) (end (aligned-length-of interval)) token-case)
  (with-slots (aligned)
      interval
    (to-string aligned :start start :end end :token-case token-case)))

(defmethod ninvert-complement :before ((interval na-alignment-interval))
  (error 'bio-sequence-op-error
         :text "Alignment intervals may not be destructively modified."))
