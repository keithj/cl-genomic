
(in-package :bio-sequence)


;;; bio-sequence generics

(defgeneric copy-sequence (bio-sequence)
  (:documentation "Returns a copy of BIO-SEQUENCE."))

(defgeneric length-of (bio-sequence)
  (:documentation "Returns the length of BIO-SEQUENCE."))

(defgeneric residue-of (bio-sequence index)
  (:documentation "Returns the residue at INDEX of BIO-SEQUENCE."))

(defgeneric (setf residue-of) (value bio-sequence index)
  (:documentation "Sets the residue at INDEX of BIO-SEQUENCE to
VALUE."))

(defgeneric to-string (bio-sequence &optional start length)
  (:documentation "Returns the string representing BIO-SEQUENCE,
starting at the first residue, or index START, and continuing to the
last residue, or for LENGTH residues."))


;;; bio-sequence io generics

(defgeneric read-fasta (input &optional callback callback-args)
  (:documentation "Reads a Fasta record from INPUT, optionally
applying function CALLBACK to the result."))

(defgeneric read-fastq (input &optional callback callback-args)
  (:documentation "Reads a Fastq record from INPUT, optionally
applying function CALLBACK to the result."))


;;; bio-graph generics

(defgeneric lookup-vertex (identity graph)
  (:documentation "Returns a vertex from GRAPH given its unique
IDENTITY."))

(defgeneric contains-vertex-p (vertex graph)
  (:documentation "Returns T if VERTEX is present in GRAPH, or NIL
otherwise."))

(defgeneric out-edges-of (vertex graph)
  (:documentation "Returns a list of the out edges of VERTEX in
GRAPH."))

(defgeneric in-edges-of (vertex graph)
  (:documentation "Returns a list of the in edges of VERTEX in
GRAPH."))

(defgeneric predecessors-of (vertex graph)
  (:documentation "Returns a list of the predecessor vertices of
VERTEX in GRAPH."))

(defgeneric successors-of (vertex graph)
  (:documentation "Returns a list of the successor vertices of VERTEX
in GRAPH."))

(defgeneric add-vertex (vertex graph)
  (:documentation "Adds VERTEX to GRAPH. Throws an error if VERTEX is
already present."))

(defgeneric remove-vertex (vertex graph)
  (:documentation "Removes VERTEX from GRAPH and removes any edges
connected to VERTEX. Throws an error if VERTEX is not present."))

(defgeneric contains-edge-p (edge graph)
  (:documentation "Returns T if EDGE is present in GRAPH, or NIL
otherwise."))

(defgeneric add-edge (source target graph)
  (:documentation "Adds an edge to GRAPH directed from vertex SOURCE
to vertex TARGET. Throws an error if an edge from SOURCE to TARGET
already exists."))

(defgeneric remove-edge (edge graph)
  (:documentation "Removes EDGE from GRAPH. Throws an error if EDGE is
not present in GRAPH."))

