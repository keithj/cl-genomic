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


;;; General purpose graph classes
(defclass graph ()
  ((vertex-table :initform (make-hash-table)
                 :reader vertex-table-of
                 :documentation "A table of all vertices in the graph,
keyed by their identity slot value.")
   (edge-table :initform (make-hash-table)
               :reader edge-table-of
               :documentation "A table of all edges in the graph,
keyed by their identity slot value."))
  (:documentation "A graph consisting of a set of indexed vertices."))

(defclass directed-graph (graph)
  ((source-table :initform (make-hash-table)
                 :reader source-table-of
                 :documentation "An table of adjacency lists keyed on
the identity slot value of source vertices. Each list contains
identities of target vertices. The table is used to follow paths in
the graph.")
   (target-table :initform (make-hash-table)
                 :reader target-table-of
                 :documentation "An table of adjacency lists keyed on
their identity slot value of target vertices. Each list contains
identities of source vertices. This table is used to look backwards
along paths in the graph.")
   (root-vertices :initform (make-hash-table)
                  :reader root-table-of
                  :documentation "A table containing the root vertices
of the graph, keyed by their identity slot value."))
  (:documentation "A directed graph consisting of set of indexed
vertices and directed edges."))

(defclass directed-acyclic-graph (directed-graph)
  ()
  (:documentation "A directed acyclic graph."))


;;; Vertex and edge classes and special behaviour methods
(defclass vertex ()
  ((identity :initarg :identity
             :reader identity-of
             :documentation "A unique fixnum identifying a vertex with
particlular semantics."))
  (:documentation "A graph vertex."))

(defclass edge ()
  ((identity :reader identity-of
             :documentation "A list of vertex identities unique to an
edge with particular semantics.")
   (source :initarg :source
           :reader source-of
           :documentation "The source vertex of the edge.")
   (target :initarg :target
           :reader target-of
           :documentation "The target vertex of the edge."))
  (:documentation "A directed graph edge from a source vertex to a
target vertex."))

;; Lazily create and cache an identity for an edge
(defmethod slot-unbound (class (obj edge) (slot (eql 'identity)))
  (with-slots (source target identity) obj
    (setf identity (list (identity-of source)
                         (identity-of target)))))


;;; Vertex methods
(defmethod add-vertex ((vertex vertex) (graph graph))
  (when (contains-vertex-p vertex graph)
    (error "Vertex ~a is already a member of graph ~a." vertex graph))
  (setf (gethash (identity-of vertex) (vertex-table-of graph)) vertex))

(defmethod add-vertex :after ((vertex vertex) (graph directed-graph))
  (setf (gethash (identity-of vertex) (root-table-of graph)) vertex))

(defmethod remove-vertex ((vertex vertex) (graph directed-graph))
  (unless (contains-vertex-p vertex graph)
    (error "Vertex ~a is not a member of graph ~a." vertex graph))
  (dolist (out-edge (out-edges-of vertex graph))
    (remove-edge out-edge graph))
  (dolist (in-edge (in-edges-of vertex graph))
    (remove-edge in-edge graph))
  (remhash (identity-of vertex) (vertex-table-of graph))
  (remhash (identity-of vertex) (root-table-of graph))
  vertex)

(defmethod remove-vertex ((identity fixnum) (graph directed-graph))
  (unless (contains-vertex-p identity graph)
    (error "Vertex identified by ~a is not a member of graph ~a."
           identity graph))
  (remove-vertex (lookup-vertex identity graph) graph))

(defmethod lookup-vertex ((identity fixnum) (graph graph))
  (gethash identity (vertex-table-of graph)))

(defmethod contains-vertex-p ((vertex vertex) (graph graph))
  (gethash (identity-of vertex) (vertex-table-of graph)))

(defmethod contains-vertex-p ((identity fixnum) (graph graph))
  (gethash identity (vertex-table-of graph)))

(defmethod out-edges-of ((source vertex) (graph directed-graph)
                         &key filter-fn)
  (let ((source-table (source-table-of graph))
        (edge-table (edge-table-of graph))
        (source-key (identity-of source)))
    (let ((edges (mapcar #'(lambda (target-key)
                             (gethash (list source-key target-key)
                                      edge-table))
                         (gethash source-key source-table))))
      (if filter-fn
          (delete-if filter-fn edges)
        edges))))

(defmethod out-edges-of ((identity fixnum) (graph directed-graph)
                         &key filter-fn)
  (out-edges-of (lookup-vertex identity graph) graph :filter-fn filter-fn))

(defmethod in-edges-of ((target vertex) (graph directed-graph)
                        &key filter-fn)
  (let ((target-table (target-table-of graph))
        (edge-table (edge-table-of graph))
        (target-key (identity-of target)))
    (let ((edges (mapcar #'(lambda (source-key)
                             (gethash (list source-key target-key)
                                      edge-table))
                         (gethash target-key target-table))))
      (if filter-fn
          (delete-if filter-fn edges)
        edges))))

(defmethod in-edges-of ((identity fixnum) (graph directed-graph)
                        &key filter-fn)
  (in-edges-of (lookup-vertex identity graph) graph :filter-fn filter-fn))

(defmethod predecessors-of ((vertex vertex) (graph directed-graph)
                            &key filter-fn)
  (let ((predecessors (collect-hash-values
                       (gethash (identity-of vertex)
                                (target-table-of graph))
                       (vertex-table-of graph))))
    (if filter-fn
        (delete-if filter-fn predecessors)
      predecessors)))

(defmethod predecessors-of ((identity fixnum) (graph directed-graph)
                            &key filter-fn)
  (predecessors-of (lookup-vertex identity graph) graph :filter-fn filter-fn))

(defmethod successors-of ((vertex vertex) (graph directed-graph)
                          &key filter-fn)
  (let ((successors (collect-hash-values
                     (gethash (identity-of vertex)
                              (source-table-of graph))
                     (vertex-table-of graph))))
    (if filter-fn
        (delete-if filter-fn successors)
      successors)))

(defmethod successors-of ((identity fixnum) (graph directed-graph)
                          &key filter-fn)
  (successors-of (lookup-vertex identity graph) graph :filter-fn filter-fn))

(defmethod root-vertices-of ((graph directed-graph))
  (loop for vertex being the hash-values in (root-table-of graph)
     collect vertex))


;;; Edge methods
(defmethod add-edge ((source vertex) (target vertex)
                     (graph directed-graph))
  (unless (contains-vertex-p source graph)
    (error "Source vertex ~a is not a member of graph ~a." source graph))
  (unless (contains-vertex-p target graph)
    (error "Target vertex ~a is not a member of graph ~a." target graph))
  (let ((source-key (identity-of source))
        (target-key (identity-of target)))
    (when (equal source-key target-key)
      (error "Cannot make an edge from vertex ~a to itself." source))
    (let ((source-table (source-table-of graph))
          (target-table (target-table-of graph))
          (new-edge (make-instance 'edge :source source :target target)))
      (when (contains-edge-p new-edge graph)
        (error (msg "An edge between vertex ~a and vertex ~a"
                    "is already present in graph ~a.") source target graph))
      (let ((out-edges (gethash source-key source-table))
            (in-edges (gethash target-key target-table))
            (new-key (identity-of new-edge)))
        (setf (gethash source-key source-table) (cons target-key out-edges)
              (gethash target-key target-table) (cons source-key in-edges)
              (gethash new-key (edge-table-of graph)) new-edge)
        (remhash target-key (root-table-of graph)))
      new-edge)))

(defmethod add-edge ((source-indentity fixnum)
                     (target-indentity fixnum) (graph directed-graph))
  (add-edge (lookup-vertex source-indentity graph)
            (lookup-vertex target-indentity graph) graph))

(defmethod add-edge ((source vertex) (target vertex)
                     (graph directed-acyclic-graph))
  (when (would-cycle-p source target graph)
    (error (msg "Adding and edge between vertices ~a and ~a"
                "would introduce a cycle in graph ~a.")
           source target graph))
  (call-next-method))

(defmethod remove-edge ((edge edge) (graph directed-graph))
  (let ((source-table (source-table-of graph))
        (target-table (target-table-of graph))
        (edge-table (edge-table-of graph))
        (source-key (identity-of (source-of edge)))
        (target-key (identity-of (target-of edge))))
    (let ((out-edge-keys (gethash source-key source-table))
          (in-edge-keys (gethash target-key target-table)))
      (unless (contains-edge-p edge graph)
        (error "Edge ~a is not present in graph ~a." edge graph))
      (update-edges source-key source-table
                    target-key target-table
                    (identity-of edge) edge-table
                    in-edge-keys out-edge-keys))
    ;; If target vertex of this edge is now absent from the
    ;; target-table this means that it has had its in-degree reduced
    ;; to 0 by the removal of the edge. This vertex is now a root.
    (unless (gethash target-key target-table)
      (setf (gethash target-key (root-table-of graph)) (target-of edge))))
  edge)

(defmethod remove-edge ((identity list) (graph directed-graph))
  (unless (contains-edge-p identity graph)
    (error "Edge identified by ~a is not a member of graph ~a."
           identity graph))
  (remove-edge (lookup-edge identity graph) graph))

(defmethod lookup-edge ((identity list) (graph graph))
  (gethash identity (edge-table-of graph)))

(defmethod contains-edge-p ((edge edge) (graph directed-graph))
  (gethash (identity-of edge) (edge-table-of graph)))

(defmethod contains-edge-p ((identity list) (graph directed-graph))
  (gethash identity (edge-table-of graph)))


;;; Graph search methods
(defmethod graph-search ((start vertex) (end vertex)
                         (graph directed-graph) &key method (test #'eql))
  (unless (contains-vertex-p start graph)
    (error "Start vertex ~a is not a member of graph ~a." start graph))
  (unless (contains-vertex-p end graph)
    (error "End vertex ~a is not a member of graph ~a." end graph))
  (let ((enqueue-fn (ecase method
                      (:depth-first #'enqueue-first)
                      (:breadth-first #'enqueue-last))))
    (graph-search-aux start end graph test
                      (list (list start))
                      enqueue-fn)))

(defmethod ancestors-of ((vertex vertex) (graph directed-acyclic-graph)
                         &key filter-fn)
  (unless (contains-vertex-p vertex graph)
    (error "Vertex ~a is not a member of graph ~a." vertex graph))
  (delete-duplicates
   (graph-traverse-aux (predecessors-of vertex graph
                                        :filter-fn filter-fn) graph
                                        #'predecessors-of filter-fn)))

(defmethod descendants-of ((vertex vertex) (graph directed-acyclic-graph)
                           &key filter-fn)
 (unless (contains-vertex-p vertex graph)
   (error "Vertex ~a is not a member of graph ~a." vertex graph))
  (delete-duplicates
   (graph-traverse-aux (successors-of vertex graph
                                      :filter-fn filter-fn) graph
                                      #'successors-of filter-fn)))

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


;;; Printing methods
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

(defmethod print-object ((edge edge) stream)
  (with-slots (source target) edge
    (format stream "<EDGE from ~a to ~a>" source target)))


;;; Internal functions
(defun roots-with-edges (directed-graph)
  "Returns a list containing all root nodes that have out-edges."
  (let ((source-table (source-table-of directed-graph)))
    (remove-if #'(lambda (vertex)
                   (not (gethash (identity-of vertex) source-table)))
               (root-vertices-of directed-graph))))

(defun would-cycle-p (source target directed-graph)
  "Returns T if adding an edge from vertex SOURCE to vertex TARGET in
DIRECTED-GRAPH would introduce a cycle."
  (some #'(lambda (root)
            (member target (graph-search root source directed-graph
                                         :method :depth-first)))
        (roots-with-edges directed-graph)))

(defun update-edges (source-key source-table
                     target-key target-table
                     edge-key edge-table in-edges out-edges)
  "Updates the adjacency list tables when EDGE is deleted. If EDGE is
the last remaining entry in the adjacency lists for vertices SOURCE
and TARGET, those records are entirely removed from the tables, rather
than leaving a hash key pointing to an empty list."
  (let ((remaining-out-edges (delete edge-key out-edges))
        (remaining-in-edges (delete edge-key in-edges)))
    (if remaining-out-edges
        (setf (gethash source-key source-table) remaining-out-edges)
      (remhash source-key source-table))
    (if remaining-in-edges
        (setf (gethash target-key target-table) remaining-out-edges)
      (remhash target-key target-table)))
  (remhash edge-key edge-table))

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

(defun collect-hash-values (keys hash-table)
  (loop for key in keys
     collect (gethash key hash-table)))
