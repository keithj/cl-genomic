
(in-package :bio-sequence)

(defgeneric copy-sequence (bio-sequence)
  (:documentation "Returns a copy of BIO-SEQUENCE."))

(defgeneric length-of (bio-sequence)
  (:documentation "Returns the length of BIO-SEQUENCE."))

(defgeneric residue-of (bio-sequence subscript)
  (:documentation "Returns the residue at index SUBSCRIPT of
BIO-SEQUENCE."))

(defgeneric (setf residue-of) (value bio-sequence subscript)
  (:documentation "Sets the residue at index SUBSCRIPT of BIO-SEQUENCE
to VALUE."))

(defgeneric to-string (bio-sequence &optional start length)
  (:documentation "Returns the string representing BIO-SEQUENCE,
starting at the first residue, or index START, and continuing to the
last residue, or for LENGTH residues."))

(defgeneric from-string (bio-sequence string)
  (:documentation "."))
