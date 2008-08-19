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
           :reader domain-of
           :documentation "The slot domain facet. A non-NIL value
indicates the type of frame of which the slot may be a member.")
   (range :initform nil
          :initarg :range
          :reader range-of
          :documentation "The slot range facet. A non-NIL value
indicates the type of values the slot may contain.")
   (value :initform nil
          :initarg :value
          :accessor value-of
          :documentation "The slot value facet."))
  (:documentation "A slot in a knowledgebase frame."))

(defclass single-valued-slot (slot)
  ()
  (:documentation "A slot whose value is a single object."))

(defclass set-valued-slot (slot)
  ()
  (:documentation "A slot whose value is a set of objects."))

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

(defclass single-valued-inverse-slot (single-valued-slot inverse-mixin)
  ())

(defclass set-valued-inverse-slot (set-valued-slot inverse-mixin)
  ())

(defclass part-of (set-valued-inverse-slot transitive-mixin)
  ((name :initform "part-of")
   (inverse :initform 'has-part)
   (inverse-name :initform "has-part")))

(defclass has-part (set-valued-inverse-slot transitive-mixin)
  ((name :initform "has-part")
   (inverse :initform 'part-of)
   (inverse-name :initform "part-of")))

(defclass subsequence-of (part-of)
  ((inverse :initform 'has-subsequence)
   (position :initarg :position
             :accessor position-of
             :documentation "The position of the subsequence in a
larger sequence.")))

(defclass has-subsequence (has-part)
  ((inverse :initform 'subsequence-of)))

(defclass instance (single-valued-slot)
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
             (format stream "Knowledgebase error~@[: ~a~]"
                     (text-of condition)))))


;; Want a shorthand for defining and adding a frame with slots


(defun framep (object)
  (subtypep (type-of object) 'frame))

(defun slotp (object)
  (subtypep (type-of object) 'slot))

(defmethod initialize-instance :after ((kb knowledgebase) &key)
  )

;;; Set default inverse-name, if possible
(defmethod initialize-instance :after ((slot inverse-mixin) &key)
  (when (and (slot-boundp slot 'inverse)
             (not (slot-boundp slot 'inverse-name)))
    (with-slots (inverse inverse-name) slot
      (setf inverse-name (string-downcase (symbol-name inverse))))))


(defgeneric contains-frame-p (frame-name &optional knowledgebase))

(defgeneric find-frame (frame-name &optional knowledgebase))
(defgeneric add-frame (frame &optional knowledgebase))
(defgeneric remove-frame (frame-name &optional knowledgebase))

(defgeneric contains-slot-p (frame slot-name))
(defgeneric find-slot (frame slot-name))
(defgeneric add-slot (frame slot))
(defgeneric remove-slot (frame slot))
(defgeneric slot-value-of (frame slot-name))
(defgeneric (setf slot-value-of) (value frame slot-name))

;;; Print a frame
(defmethod print-object ((frame frame) stream)
  (let ((*print-circle* t))
    (with-slots (name slots) frame
      (format stream "<FRAME ~a ~a>" name slots))))

(defmethod print-object ((slot slot) stream)
  (let ((*print-circle* t))
    (with-slots (name value) slot
      (format stream "<SLOT ~a: ~a>" name value))))

;;; Frame methods
(defmethod find-frame ((frame-name string)
                       &optional (kb *default-knowledgebase*))
  (multiple-value-bind (frame present-p)
      (gethash frame-name (frames-of kb))
    (unless present-p
      (error 'knowledgebase-error :text
             (format nil "~a is not present in ~a." frame-name kb)))
    frame))

(defmethod contains-frame-p ((frame-name string)
                             &optional (kb *default-knowledgebase*))
  (gethash frame-name (frames-of kb)))

(defmethod contains-frame-p ((frame frame) &optional
                             (kb *default-knowledgebase*))
  (gethash (name-of frame) (frames-of kb)))

(defmethod add-frame :before ((frame frame) &optional
                              (kb *default-knowledgebase*))
  (when (contains-frame-p frame kb)
    (error 'knowledgebase-error :text
           (format nil "~a is already present in ~a." frame kb))))

(defmethod add-frame ((frame frame) &optional
                      (kb *default-knowledgebase*))
  (setf (gethash (name-of frame) (frames-of kb)) frame))

(defmethod remove-frame :before (frame &optional
                                 (kb *default-knowledgebase*))
  (unless (contains-frame-p frame kb)
    (error 'knowledgebase-error :text
           (format nil "~a is not present in ~a." frame kb))))

(defmethod remove-frame ((frame frame)
                         &optional (kb *default-knowledgebase*))
  (remhash (name-of frame) (frames-of kb))
  frame)


;;; Slot methods
(defmethod find-slot :before ((frame frame) (slot-name string))
  (unless (contains-slot-p frame slot-name)
    (error 'knowledgebase-error :text
           (format nil "~a is not a slot of ~a." slot-name frame))))

(defmethod find-slot ((frame frame) (slot-name string))
  (let* ((slots (slots-of frame))
         (slot-pos (position slot-name slots :key #'name-of :test #'equal)))
    (aref slots slot-pos)))

(defmethod contains-slot-p ((frame frame) (slot-name string))
  (find slot-name (slots-of frame) :key #'name-of :test #'equal))

(defmethod contains-slot-p ((frame frame) (slot slot))
  (find slot (slots-of frame)))

(defmethod add-slot :before ((frame frame) (slot slot))
  (when (contains-slot-p frame slot)
    (error 'knowledgebase-error :text
           (format nil "~a is already a slot of ~a." slot frame)))
  (let ((slot-domain (domain-of slot)))
    (when slot-domain
      (unless (subtypep (type-of frame) slot-domain)
        (error 'knowledgebase-error :text
               (format nil (msg "Invalid assignment of ~a into ~a:"
                                "the domain of ~a is restricted"
                                "to subtypes of ~a.")
                       slot frame slot slot-domain))))))

(defmethod add-slot ((frame frame) (slot slot))
  (vector-push-extend slot (slots-of frame))
  frame)

(defmethod remove-slot :before ((frame frame) (slot slot))
  (unless (contains-slot-p frame slot)
    (error 'knowledgebase-error :text
           (format nil "~a is not a slot of ~a." slot frame))))

(defmethod remove-slot ((frame frame) (slot slot))
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


;;; Slot value methods
(defmethod slot-value-of ((frame frame) (slot-name string))
  (unless (contains-slot-p frame slot-name)
    (error 'knowledgebase-error :text
           (format nil "~a is not a slot of ~a." slot-name frame)))
  ;; Need to delegate correctly so that dispatch on slot type can
  ;; occur
  (slot-value-of frame (find-slot frame slot-name)))

(defmethod slot-value-of ((frame frame) (slot single-valued-slot))
  (value-of slot))

;; Need to separate set-valued and single-valued
;; (defmethod slot-value-of ((frame frame) (slot transitive-mixin))
;;   (let ((value (value-of slot)))
;;     (cond ((null value)
;;            nil)
;;           ((framep value)
;;            (cons value (slot-value-of value (name-of slot))))
;;           ((listp value)
;;            (loop
;;               for subval in value
;;               nconc (if (framep subval)
;;                         (append (list subval) (slot-value)))
;;                 ))))
;;   )

(defmethod (setf slot-value-of) (value (frame frame) (slot-name string))
  (unless (contains-slot-p frame slot-name)
    (error 'knowledgebase-error :text
           (format nil "~a is not a slot of ~a." slot-name frame)))
  (update-slot-value frame (find-slot frame slot-name) value))



;; Possibly separate methods for direct slot values transitive slot
;; values?

;;; Slot update methods
(defmethod update-slot-value :before ((frame frame) (slot single-valued-slot)
                                      value)
  (let ((slot-range (range-of slot)))
    (when slot-range
      (unless (subtypep (type-of value) slot-range)
        (error 'knowledgebase-error :text
               (format nil (msg "Invalid slot value ~a: the range of"
                                "~a is restricted to subtypes of ~a.")
                       value slot slot-range))))))

(defmethod update-slot-value ((frame frame) (slot single-valued-slot)
                              value)
  (setf (value-of slot) value))

(defmethod update-slot-value :before ((frame frame) (slot set-valued-slot)
                                      (values list))
  (let ((slot-range (range-of slot)))
    (when slot-range
      (dolist (value values)
        (unless (subtypep (type-of value) slot-range)
          (error 'knowledgebase-error :text
                 (format nil (msg "Invalid slot value ~a: the range of"
                                  "~a is restricted to subtypes of ~a.")
                         value slot slot-range)))))))

(defmethod update-slot-value ((frame frame) (slot set-valued-slot)
                              (values list))
  (setf (value-of slot) values))

(defmethod update-slot-value ((subject frame) (slot single-valued-inverse-slot)
                              (object frame))
  (let* ((inverse-name (inverse-name-of slot))
         (inverse-slot (find-or-make-slot object inverse-name
                                          (inverse-of slot)
                                          :name inverse-name
                                          :domain (range-of slot)
                                          :range (domain-of slot))))
    (update-inverse-slot-value subject inverse-slot object))
  (call-next-method))

(defmethod update-slot-value ((subject frame) (slot set-valued-inverse-slot)
                              (objects list))
  (let ((inverse-name (inverse-name-of slot)))
    (dolist (object objects)
      (let ((inverse-slot (find-or-make-slot object inverse-name
                                             (inverse-of slot)
                                             :name inverse-name
                                             :domain (range-of slot)
                                             :range (domain-of slot))))
        (update-inverse-slot-value subject inverse-slot object))))
  (call-next-method))

(defmethod update-inverse-slot-value ((object frame) (slot single-valued-slot)
                                      (subject frame))
  (setf (value-of slot) object))

(defmethod update-inverse-slot-value ((object frame) (slot set-valued-slot)
                                      (subject frame))
  ;; FIXME -- should set-valued slots be using add-slot-value rather
  ;; than setting the whole set?
  (let ((values (value-of slot)))
    (setf (value-of slot) (adjoin object values))))


(defun find-or-make-slot (frame slot-name slot-class &rest slot-initargs)
  (if (contains-slot-p frame slot-name)
      (find-slot frame slot-name)
    (let ((slot (apply #'make-instance slot-class slot-initargs)))
      (add-slot frame slot)
      slot)))

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
