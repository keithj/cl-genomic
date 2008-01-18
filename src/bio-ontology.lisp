
(in-package :bio-sequence)


;;; Ontology term (vertex) and relationship (edge) classes and special
;;; behaviour methods

(defclass ontology-term (vertex)
  ((name :initform nil
         :initarg :name
         :reader name-of
         :documentation "The human-readable name of the term.")
   (definition :initarg :definition
               :reader definition-of
               :documentation "The text defining what the term
represents." )
   (synonyms :initform nil
             :initarg :synonyms
             :accessor synonyms-of
             :documentation "A list of names synonymous with term
name.")
   (comment :initform nil
            :initarg :comment
            :accessor comment-of
            :documentation "A text comment which does not affect the
semantics of the term.")))

(defclass attributed-mixin ()
  ((attributes :initform nil
               :initarg :attributes
               :accessor attributes-of
               :documentation "An alist of attributes with keywords as
keys."))
  (:documentation "A mixin for classes with a small number of
arbitrary attributes which are not formally adopted as slots."))

(defclass ontology-instance-mixin (ontology-term attributed-mixin)
  ((ontology-class :initarg :ontology-class
                   :reader ontology-class-of
                   :documentation "The identity of an ontology term
which defines a class of which this vertex represents an instance."))
  (:documentation "A graph vertex which represents the instantiation
of an ontological class."))

(defclass instance-relationship (edge attributed-mixin)
  ((predicate :initarg :predicate
              :reader predicate-of))
  (:documentation "A graph edge describing a relationship between two
ontology instances."))


;; Alias source to subject and target to object as a more appropriate
;; name for an ontology
(defmethod subject-of ((obj instance-relationship))
  (slot-value obj 'source))

(defmethod (setf subject-of) ((obj instance-relationship)
                              (term ontology-term))
  (setf (slot-value obj 'source) term))

(defmethod object-of ((obj instance-relationship))
  (slot-value obj 'target))

(defmethod (setf object-of) ((obj instance-relationship)
                             (term ontology-term))
  (setf (slot-value obj 'target) term))

;; Lazily create and cache an identity for a relationship
(defmethod slot-unbound (class (obj instance-relationship)
                         (slot (eql 'identity)))
  (with-slots (source predicate target identity) obj
    (setf identity (list (identity-of source)
                         (identity-of predicate)
                         (identity-of target)))))

(defmethod print-object ((term ontology-term) stream)
  (format stream "<ONTOLOGY-TERM ~a>" (slot-value term 'identity)))

(defmethod print-object ((relationship instance-relationship) stream)
  (with-slots (source predicate target) relationship
    (format stream "<ONTOLOGY-RELATIONSHIP ~a ~a ~a>"
            source predicate target)))

(defmethod attribute-of ((obj attributed-mixin) key)
  (with-slots (attributes) obj
    (assocdr key attributes)))

(defmethod (setf attribute-of) (value (obj attributed-mixin) key)
  (with-slots (attributes) obj
    (if (assoc key attributes)
        (rplacd (assoc key attributes) value)
      (setf attributes (acons key value attributes))))
  value)
