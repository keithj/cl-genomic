;;;
;;; Copyright (C) 2008, Keith James. All rights reserved.
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

;; To do:

;; Limit the methods to accepting frame or slot objects rather than
;; their names - saves having to write many methods twice. If all we
;; have are the names, just look up the frame or slot by name first.

;; Inverse slots need special value setting methods to set up the
;; inverse relationship

;; Slots with domain or range specified need methods to check the
;; validity of those constraints (the OO way of looking at it).

;; Transitive slots need methods for collecting transitive closure.

;; Reflexive slots need methods for adding the frame to the values.

(defparameter *default-knowledgebase* nil)

(defclass knowledgebase ()
  ((frames :initform (make-hash-table :test #'equal)
           :reader frames-of
           :documentation "A table of all frames in the knowledgebase,
keyed by their name."))
  (:documentation "A knowledgebase consisting of a graph of indexed
frames."))

(defclass slot ()
  ((name :initarg :name
         :reader name-of
         :documentation "The slot name.")
   (domain :initform nil
           :initarg :domain
           :accessor domain-of
           :documentation "The slot domain facet.")
   (range :initform nil
          :initarg :range
          :accessor range-of
          :documentation "The slot range facet.")
   (value :initform nil
          :initarg :value
          :accessor value-of
          :documentation "The slot value facet."))
  (:documentation "A slot in a knowledgebase frame."))

(defclass reflexive-mixin ()
  ()
   (:documentation "A mixin identifying a relationship that is
reflexive."))

(defclass transitive-mixin ()
  ()
  (:documentation "A mixin identifying a relationship that is
transitive."))

(defclass inverse-mixin ()
  ((inverse :initarg :inverse
            :reader inverse-of
            :documentation "The relationship class of which this is
the inverse.")
   (inverse-name :initarg :inverse-name
                 :reader inverse-name-of
                 :documentation "The name of the slot of which this is
the inverse."))
  (:documentation "A slot which has an inverse relationship with
another."))

(defclass part-of (slot inverse-mixin transitive-mixin)
  ((name :initform "part-of")
   (inverse :initform 'has-part)))

(defclass has-part (slot inverse-mixin transitive-mixin)
  ((name :initform "has-part")
   (inverse :initform 'part-of)))

(defclass subsequence-of (part-of)
  ((inverse :initform 'has-subsequence)
   (position :initarg :position
             :accessor position-of
             :documentation "The position of the subsequence in a
larger sequence.")))

(defclass has-subsequence (has-part)
  ((inverse :initform 'subsequence-of)))

(defclass instance (slot)
  ((ontology-class :initarg :onto-class
                   :reader onto-class-of
                   :documentation "The ontology class of which the
slot indicates instantiation."))
  (:documentation "A slot denoting instantiation of an ontology
class."))

(defclass frame ()
  ((name :initarg :name
         :reader name-of
         :documentation "A name, unique within the namespace of a
knowledgebase, identifying a frame with particlular semantics.")
   (slots :initform (make-array 5 :element-type 'slot :adjustable t
                                :fill-pointer 0)
          :initarg :slots
          :accessor slots-of
          :documentation "A sequence of slots."))
  (:documentation "A knowledgebase frame. Since slot access requires
linear time, it is best not to have more than about 10 slots per
frame."))

(define-condition knowledgebase-error (error)
  ((text :initform nil
         :initarg :text
         :reader text-of
         :documentation "Error message text."))
  (:report (lambda (condition stream)
             (format stream "Knowledgebase error~@[: ~a~]."
                     (text-of condition)))))


;; Want a shorthand for defining and adding a frame with slots


(defmethod initialize-instance :after ((kb knowledgebase) &key)
  )

(defmethod initialize-instance :after ((slot inverse-mixin) &key)
  (with-slots (inverse inverse-name) slot
    (setf inverse-name (string-downcase (symbol-name inverse)))))


;;; Print a frame
(defmethod print-object ((frame frame) stream)
  (with-slots (name slots) frame
      (format stream "<FRAME ~a ~a>" name slots)))

(defmethod print-object ((slot slot) stream)
  (with-slots (name) slot
    (format stream "<SLOT ~a>" name)))


;;; Frame methods with string parameters
(defmethod contains-frame-p ((name string)
                             &optional (kb *default-knowledgebase*))
  (gethash name (frames-of kb)))

(defmethod find-frame ((name string)
                       &optional (kb *default-knowledgebase*))
  (multiple-value-bind (frame present-p)
      (gethash name (frames-of kb))
    (unless present-p
      (error 'knowledgebase-error :text
             (format nil "~a is not present in ~a." name kb)))
    frame))

(defmethod remove-frame ((name string)
                         &optional (kb *default-knowledgebase*))
  (remove-frame (find-frame name kb)))

(defmethod contains-slot-p ((frame frame) (name string)
                            &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (find name (slots-of frame) :key #'name-of :test #'equal))


;;; Slot method with string parameters
(defmethod find-slot :before ((frame frame) (name string)
                              &optional (kb *default-knowledgebase*))
  (unless (contains-slot-p frame name)
    (error 'knowledgebase-error :text
           (format nil "~a is not a slot of ~a in ~a." name frame kb))))

(defmethod find-slot ((frame frame) (name string)
                      &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (let* ((slots (slots-of frame))
         (slot-pos (position name slots :key #'name-of :test #'equal)))
    (aref slots slot-pos)))

(defmethod slot-value-of ((frame-name string) (slot-name string)
                          &optional (kb *default-knowledgebase*))
  (value-of (find-slot (find-frame frame-name kb) slot-name kb)))

(defmethod (setf slot-value-of) (value (frame-name string) (slot-name string)
                                 &optional (kb *default-knowledgebase*))
  (setf (value-of (find-slot (find-frame frame-name kb) slot-name kb))
        value))


;;; Frame methods
(defmethod contains-frame-p ((frame frame)
                             &optional (kb *default-knowledgebase*))
  (gethash (name-of frame) (frames-of kb)))

(defmethod add-frame :before ((frame frame)
                              &optional (kb *default-knowledgebase*))
  (when (contains-frame-p frame kb)
    (error 'knowledgebase-error :text
           (format nil "~a is already present in ~a." frame kb))))

(defmethod add-frame ((frame frame)
                      &optional (kb *default-knowledgebase*))
  (setf (gethash (name-of frame) (frames-of kb)) frame))

(defmethod remove-frame :before (frame
                                 &optional (kb *default-knowledgebase*))
  (unless (contains-frame-p frame kb)
    (error 'knowledgebase-error :text
           (format nil "~a is not present in ~a." frame kb))))

(defmethod remove-frame ((frame frame)
                         &optional (kb *default-knowledgebase*))
  (remhash (name-of frame) (frames-of kb))
  frame)


;;; Slot methods
(defmethod contains-slot-p ((frame frame) (slot slot)
                            &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (find slot (slots-of frame)))

(defmethod add-slot :before ((frame frame) (slot slot)
                             &optional (kb *default-knowledgebase*))
  (unless (contains-frame-p frame kb)
    (error 'knowledgebase-error :text
           (format nil "~a is not present in ~a." frame kb)))
  (when (contains-slot-p frame slot kb)
    (error 'knowledgebase-error :text
           (format nil "~a is already a slot of ~a in ~a."
                   slot frame kb))))

(defmethod add-slot ((frame frame) (slot slot)
                     &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (vector-push-extend slot (slots-of frame))
  frame)

(defmethod remove-slot :before ((frame frame) (slot slot)
                                &optional (kb *default-knowledgebase*))
  (unless (contains-frame-p frame kb)
    (error 'knowledgebase-error :text
           (format nil "~a is not present in ~a." frame kb)))
  (unless (contains-slot-p frame slot kb)
    (error 'knowledgebase-error :text
           (format nil "~a is not a slot of ~a in ~a." slot frame kb))))

(defmethod remove-slot ((frame frame) (slot slot)
                        &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (let* ((slots (slots-of frame))
         (last-slot-pos (1- (length slots)))
         (slot-pos (position slot slots)))
    ;; We aim to delete a slot by adjusting the fill-pointer. If the
    ;; slot to be removed isn't the last one, copy the last one over
    ;; it, then remove the last one, which is now a duplicate.
    (when (\= slot-pos last-slot-pos)
      (setf (aref slots slot-pos) (aref slots last-slot-pos)))
    (adjust-array slots (length slots) :fill-pointer last-slot-pos))
  frame)

;;; Slot removal currently does not remove inverse slots - should it?


(defmethod slot-value-of ((frame frame) (slot slot)
                          &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (value-of slot))

(defmethod (setf slot-value-of) (value (frame frame) (slot slot)
                                 &optional (kb *default-knowledgebase*))
  (declare (ignore kb))
  (setf (value-of slot) value))

(defmethod (setf slot-value-of) ((object frame) (subject frame)
                                 (slot inverse-mixin)
                                 &optional (kb *default-knowledgebase*))
  (let ((inverse-name (inverse-name-of slot)))
    (if (contains-slot-p object inverse-name)
        (setf (slot-value-of object inverse-name) subject)
      (let ((inverse-slot (make-instance (inverse-of slot)
                                         :name inverse-name
                                         :domain (range-of slot)
                                         :range (domain-of slot)
                                         :value subject)))
        (add-slot object inverse-slot kb))))
  (call-next-method))

(defmethod slot-value-of ((frame frame) (slot transitive-mixin)
                          &optional (kb *default-knowledgebase*))
  ;; apply slot-value-of to all values of slot that are frames

  ;; all values of the slot should be of the same type (or if
  ;; set-valued, all members of the slot's set)

  ;; Need to check domain and range to ensure that this is true
  )



(defun graph-traverse-aux (sub-graph dag traversal-fn filter-fn)
  "Traverses directed acyclic graph DAG, starting from SUB-GRAPH (a
vertex or list of a vertex's successors or predecessors), using the
traversal function TRAVERSAL-FN to find the next successors or
predecessors. Suitable traversal functions are PREDECESSORS-OF and
SUCCESSORS-OF. The vertices traversed are returned in a list which may
include duplicates where DAG contains diamonds."
  (cond ((null sub-graph)
         nil)
        ((atom sub-graph)
         (cons sub-graph
               (graph-traverse-aux (funcall traversal-fn sub-graph
                                            dag :filter-fn filter-fn)
                                   dag traversal-fn filter-fn)))
        (t
         (concatenate 'list
                      (graph-traverse-aux (car sub-graph)
                                          dag traversal-fn filter-fn)
                      (graph-traverse-aux (cdr sub-graph)
                                          dag traversal-fn filter-fn)))))

