
(in-package :bigvector)

(defgeneric length-of (bigvector)
  (:documentation "Returns the number of elements in BIGVECTOR."))

(defgeneric element-type-of (bigvector)
  (:documentation "Returns the type of elements in BIGVECTOR."))

(defgeneric bvref (bigvector subscript)
  (:documentation "Accesses the BIGVECTOR element specified by
SUBSCRIPT."))

(defgeneric bigvector-subseq (bigvector start &optional end)
  (:documentation "Creates a bigvector that is a copy of the
subsequence of BIGVECTOR bounded by START and END. START specifies an
offset into BIGVECTOR and marks the beginning position of the
subsequence, END marks the position following the last element of the
subsequence. If END is not supplied, its value defaults to the length
of BIGVECTOR. This function is analagous to SUBSEQ for sequences."))

(defgeneric bigvector-concatenate (&rest bigvectors)
  (:documentation "Returns a new bigvector that contains all the
elements of BIGVECTORS in the order that they are supplied. The
element-type of the new bigvector is the same as that of the first
argument."))

(defgeneric bigvector-copy (source dest &key source-start source-end
                                   dest-start key)
  (:documentation "Destructively copies a subsequence from SOURCE that
is bounded by SOURCE-START and SOURCE-END into DEST. DEST-START
specifies the offset into DEST at which the subsequence is
inserted. The optional parameter KEY specifies a transformation
function to be applied to each element in the subsequence before its
insertion into DEST."))
