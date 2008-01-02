
(in-package :bio-sequence)

;; General purpose graph classes

(defclass graph ()
  ((vertex-table :initarg vertex-table
                 :initform (make-hash-table :test #'equal)
                 :reader vertex-table-of
                 :documentation "A table of all vertices in the graph,
indexed by their unique identifying string."))
  (:documentation "A graph consisting of a set of indexed vertices."))

(defclass directed-graph (graph)
  ((source-table :initform (make-hash-table)
                 :reader source-table-of
                 :documentation "An table of adjacency lists keyed on
source vertices. Each list contains edges. The table is used to follow
paths in the graph.")
   (target-table :initform (make-hash-table)
                 :reader target-table-of
                 :documentation "An table of adjacency lists keyed on
target vertices. Each list contains edges. This table is used to look
backwards along paths in the graph.")
   (root-vertices :initform nil
                  :accessor root-vertices-of
                  :documentation "A list containing the root vertices
of the graph."))
  (:documentation "A directed graph consisting of set of indexed
vertices and directed edges."))

(defclass directed-acyclic-graph (directed-graph)
  ()
  (:documentation "A directed acyclic graph."))


;; Vertex and edge classes

(defclass vertex ()
  ((identity :initarg :identity
             :reader identity-of
             :documentation "A unique string identifying the vertex."))
  (:documentation "A graph vertex."))

(defclass edge ()
  ((source :initarg :source
           :reader source-of
           :documentation "The source vertex of the edge.")
   (target :initarg :target
           :reader target-of
           :documentation "The target vertex of the edge."))
  (:documentation "A directed graph edge from a source vertex to a
target vertex."))

(defclass ontology-term (vertex)
  ((name :initarg :name
         :reader name-of)
   (definition :initarg :definition
               :reader definition-of)
   (comment :initarg :comment
            :accessor comment-of)))


;; Triple class

(defclass triple (edge)
  ((predicate :initarg :predicate
              :reader predicate-of)))

; Initializer and accessor aliases for triple slots

(defmethod initialize-instance :after ((obj triple) &key subject object)
  (when subject
    (setf (slot-value obj 'source) subject))
  (when object
    (setf (slot-value obj 'target) object)))

(defmethod subject-of ((obj triple))
  (slot-value obj 'source))

(defmethod (setf subject-of) ((obj triple) (term ontology-term))
  (setf (slot-value obj 'source) term))

(defmethod object-of ((obj triple))
  (slot-value obj 'target))

(defmethod (setf object-of) ((obj triple) (term ontology-term))
  (setf (slot-value obj 'target) term))


;; Vertex methods

(defmethod add-vertex ((vertex vertex) (graph graph))
  (let ((vertex-table (vertex-table-of graph))
        (identity (identity-of vertex)))
    (if (gethash identity vertex-table)
        (error "vertex ~a is already a member of graph ~a" vertex graph)
      (setf (gethash identity vertex-table) vertex))))

(defmethod add-vertex :after ((vertex vertex) (graph directed-graph))
  (push vertex (root-vertices-of graph)))

(defmethod remove-vertex ((vertex vertex) (graph directed-graph))
  (unless (contains-vertex-p vertex graph)
    (error "from vertex ~a is not a member of graph ~a" vertex graph))
  (dolist (out-edge (out-edges-of vertex graph))
    (remove-edge out-edge graph))
  (dolist (in-edge (in-edges-of vertex graph))
    (remove-edge in-edge graph))
  (remhash (identity-of vertex) (vertex-table-of graph))
  (setf (root-vertices-of graph)
        (delete vertex (root-vertices-of graph)))
  vertex)

(defmethod lookup-vertex ((identity string) (graph graph))
  (gethash identity (vertex-table-of graph)))

(defmethod contains-vertex-p ((vertex vertex) (graph graph))
  (gethash (identity-of vertex) (vertex-table-of graph)))

(defmethod contains-vertex-p ((identity string) (graph graph))
  (gethash identity (vertex-table-of graph)))

(defmethod out-edges-of ((source vertex) (graph directed-graph))
  (gethash source (source-table-of graph)))

(defmethod out-edges-of ((identity string) (graph directed-graph))
  (gethash (lookup-vertex identity graph) (source-table-of graph)))

(defmethod in-edges-of ((target vertex) (graph directed-graph))
  (gethash target (target-table-of graph)))

(defmethod in-edges-of ((identity string) (graph directed-graph))
  (gethash (lookup-vertex identity graph) (target-table-of graph)))

(defmethod predecessors-of ((vertex vertex) (graph directed-graph))
  (mapcar #'source-of (in-edges-of vertex graph)))

(defmethod predecessors-of ((identity string) (graph directed-graph))
  (mapcar #'source-of (in-edges-of identity graph)))

(defmethod successors-of ((vertex vertex) (graph directed-graph))
  (mapcar #'target-of (out-edges-of vertex graph)))

(defmethod successors-of ((identity string) (graph directed-graph))
  (mapcar #'target-of (out-edges-of identity graph)))


;; Edge methods

(defmethod add-edge ((source vertex) (target vertex)
                     (graph directed-graph))
  (unless (contains-vertex-p source graph)
    (error "from vertex ~a is not a member of graph ~a" source graph))
  (unless (contains-vertex-p target graph)
    (error "to vertex ~a is not a member of graph ~a" target graph))
  (when (string= (identity-of source) (identity-of target))
    (error "cannot make an edge from vertex ~a to itself" source))
  (let ((source-table (source-table-of graph))
        (target-table (target-table-of graph)))
    (let ((out-edges (gethash source source-table))
          (in-edges (gethash target target-table)))
      (if (find target out-edges :key #'target-of)
          (error "an edge between vertex ~a and vertex ~a is already
present in graph ~a" source target graph)
        (let ((out-edge (make-instance 'edge :source source :target target)))
          (setf (gethash source source-table) (cons out-edge out-edges)
                (gethash target target-table) (cons out-edge in-edges)
                (root-vertices-of graph)
                (delete target (root-vertices-of graph)))
          out-edge)))))

(defmethod add-edge ((source-indentity string)
                     (target-indentity string) (graph directed-graph))
  (add-edge (lookup-vertex source-indentity graph)
            (lookup-vertex target-indentity graph) graph))

(defmethod add-edge ((source vertex) (target vertex)
                     (graph directed-acyclic-graph))
  (when (would-cycle-p source target graph)
    (error "adding and edge between vertices ~a and ~a would introduce
a cycle in graph ~a" source target graph))
  (call-next-method))

(defmethod remove-edge ((edge edge) (graph directed-graph))
  (let ((source-table (source-table-of graph))
        (target-table (target-table-of graph))
        (source (source-of edge))
        (target (target-of edge)))
    (let ((out-edges (gethash source source-table))
          (in-edges (gethash target target-table)))
      (if (find edge out-edges)
          (update-edges source source-table
                        target target-table edge in-edges out-edges)
        (error "edge ~a is not present in graph ~a" edge graph)))
    ;; If target vertex of this egde is now absent from the
    ;; target-table this means that it has had its in-degree reduced
    ;; to 0 by the removal of the edge. This vertex is now a root.
    (setf (root-vertices-of graph)
          (cons target (root-vertices-of graph))))
  edge)

(defmethod remove-edge ((identities cons) (graph directed-graph))
  (let ((source-indentity (car identities))
        (target-indentity (cdr identities)))
    (unless (and (subtypep (type-of source-indentity) 'string)
                 (subtypep (type-of target-indentity) 'string))
      (error "invalid identities argument ~a; must be a dotted pair of
strings identifying a pair of vertices" identities))
    (let ((source-table (source-table-of graph))
          (target-table (target-table-of graph))
          (source (lookup-vertex source-indentity graph))
          (target (lookup-vertex target-indentity graph)))
      (let ((out-edges (gethash source source-table))
            (in-edges (gethash target target-table)))
        (let ((edge (find-if
                     #'(lambda (edge)
                         (and (string= source-indentity
                                       (identity-of (source-of edge)))
                              (string= target-indentity
                                       (identity-of (target-of edge)))))
                     out-edges)))
          (if edge
              (update-edges source source-table
                            target target-table edge in-edges out-edges)
            (error "edge between vertices ~a and ~a is not present in
graph ~a" source-indentity target-indentity graph))
          ;; If target vertex of this egde is now absent from the
          ;; target-table this means that it has had its in-degree reduced
          ;; to 0 by the removal of the edge. This vertex is now a root.
          (setf (root-vertices-of graph)
                (cons target (root-vertices-of graph)))
          edge)))))

(defmethod contains-edge-p ((edge edge) (graph directed-graph))
  (find edge (out-edges-of (source-of edge) graph)))


;; Graph search methods

(defmethod graph-search ((start vertex) (end vertex)
                         (graph directed-graph) &key method (test #'eql))
  (unless (contains-vertex-p start graph)
    (error "from vertex ~a is not a member of graph ~a" start graph))
  (unless (contains-vertex-p end graph)
    (error "to vertex ~a is not a member of graph ~a" end graph))
  (let ((enqueue-fn (ecase method
                      (:depth-first #'enqueue-first)
                      (:breadth-first #'enqueue-last))))
    (graph-search-aux start end graph test
                      (list (list start))
                      enqueue-fn)))

(defmethod ancestors ((vertex vertex) (graph directed-graph))
  (delete-duplicates
   (graph-traverse-aux (predecessors-of vertex graph) graph
                       #'predecessors-of)))

(defmethod descendants ((vertex vertex) (graph directed-graph))
  (delete-duplicates
   (graph-traverse-aux (successors-of vertex graph) graph
                       #'successors-of)))

(defun graph-traverse-aux (sub-graph directed-graph traversal-fn)
  (cond ((null sub-graph)
         nil)
        ((atom sub-graph)
         (cons sub-graph
               (graph-traverse-aux (funcall traversal-fn
                                            sub-graph directed-graph)
                                   directed-graph traversal-fn)))
        (t
         (concatenate 'list
                      (graph-traverse-aux (car sub-graph)
                                          directed-graph traversal-fn)
                      (graph-traverse-aux (cdr sub-graph)
                                          directed-graph traversal-fn)))))

;; Printing methods

(defmethod print-object ((graph graph) stream)
  (format stream "<GRAPH ~a vertices>"
          (hash-table-count (slot-value graph 'vertex-table))))

;; These graph summaries are expensive, so avoid using them for
;; logging
(defmethod print-object ((graph directed-graph) stream)
  (format stream "<DIRECTED-GRAPH ~a vertices, ~a edge~:p>"
          (hash-table-count (slot-value graph 'vertex-table))
          (loop for edges being the hash-values of
               (slot-value graph 'source-table)
             summing (length edges) into edge-count
             finally (return edge-count))))

(defmethod print-object ((graph directed-acyclic-graph) stream)
  (format stream "<DIRECTED-ACYCLIC-GRAPH ~a vertices, ~a edge~:p>"
          (hash-table-count (slot-value graph 'vertex-table))
          (loop for edges being the hash-values of
               (slot-value graph 'source-table)
             summing (length edges) into edge-count
             finally (return edge-count))))

(defmethod print-object ((vertex vertex) stream)
  (format stream "<VERTEX ~a>" (slot-value vertex 'identity)))

(defmethod print-object ((term ontology-term) stream)
  (format stream "<ONTOLOGY-TERM ~a>" (slot-value term 'identity)))

(defmethod print-object ((edge edge) stream)
  (format stream "<EDGE from ~a to ~a>"
          (slot-value edge 'source) (slot-value edge 'target)))

(defmethod print-object ((triple triple) stream)
  (format stream "<TRIPLE ~a ~a ~a>"
          (slot-value triple 'source)
          (slot-value triple 'predicate)
          (slot-value triple 'target)))


;; Internal functions

(defun roots-with-edges (directed-graph)
  "Returns a list containing all root nodes that have out-edges."
  (let ((source-table (source-table-of directed-graph)))
    (remove-if #'(lambda (vertex)
                   (not (gethash vertex source-table)))
               (root-vertices-of directed-graph))))

(defun would-cycle-p (source target directed-graph)
  "Returns T if adding an edge from vertex SOURCE to vertex TARGET in
DIRECTED-GRAPH would introduce a cycle."
  (some #'(lambda (root)
            (member target (graph-search root source directed-graph
                                         :method :depth-first)))
        (roots-with-edges directed-graph)))

(defun update-edges (source source-table
                     target target-table edge in-edges out-edges)
  "Updates the adjacency list tables when EDGE is deleted. If EDGE is
the last remaining entry in the adjacency lists for vertices SOURCE
and TARGET, those records are entirely removed from the tables, rather
than leaving a hash key pointing to an empty list."
  (let ((remaining-out-edges (delete edge out-edges))
        (remaining-in-edges (delete edge in-edges)))
    (if remaining-out-edges
        (setf (gethash source source-table) remaining-out-edges)
      (remhash source source-table))
    (if remaining-in-edges
        (setf (gethash target source-table) remaining-out-edges)
      (remhash target target-table))))

(defun graph-search-aux (start end graph test queue enqueue-fn)
  "Searches GRAPH for a path between vertices START and END, using
TEST to compare vertices for equivalence. The paths currently being
followed in search of END are held in list QUEUE. Newly extended paths
are enqueued by ENQUEUE-FN to allow different search strategies to be
parameterized."
  (cond ((endp queue)
         nil)
        ((funcall test end (first (first queue)))
         (nreverse (first queue)))
        (t
         (graph-search-aux start end graph test
                           (funcall enqueue-fn
                                    (extend-paths (first queue) graph)
                                    (rest queue)) 
                           enqueue-fn))))

(defun extend-paths (path graph)
  "Returns a list of all possible paths made by prepending the
successors of the first vertex in PATH, a list of vertices, to PATH."
  (mapcar #'(lambda (successor)
              (cons successor path))
          (remove-if #'(lambda (successor)
                         (member successor path))
                     (successors-of (first path) graph))))

(defun enqueue-first (paths queue)
  "Appends PATHS to the front of QUEUE, thereby creating the
conditions for a breadth-first search."
  (append paths queue))

(defun enqueue-last (paths queue)
  "Appends PATHS to the end of QUEUE, thereby creating the conditions
for a depth-first search."
  (append queue paths))


