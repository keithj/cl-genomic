
(in-package :bio-sequence)

(defvar +gff3-dot-column+ "."
  "The GFF3 empty column string.")
(defvar +gff3-forward-strand+ "+"
  "The GFF3 forward strand string.")
(defvar +gff3-reverse-strand+ "-"
  "The GFF3 reverse strand string.")
(defvar +gff3-unknown-strand+ "?"
  "The GFF3 unknown strand string.")

(defvar *gff3-unsafe-chars*
  '(#\; #\= #\% #\& #\,)
  "Characters that must be escaped in GFF3 attribute context.")
(defvar *gff3-valid-id-chars*
  (concatenate 'list
               (loop for i from 48 to 57 collect (code-char i))
               (loop for i from 65 to 90 collect (code-char i))
               (loop for i from 97 to 122 collect (code-char i))
               '(#\. #\: #\^ #\* #\$ #\@ #\! #\+ #\_ #\? #\- #\|))
  "Characters that are legal without escapes in GFF3 seqid context.")

(defvar *gff3-reserved-attr-tags*
  '("ID" "Name" "Alias" "Parent" "Target" "Gap" "Derives_from"
    "Note" "Dbxref" "Ontology_term" "Index")
  "Reserved GFF3 attribute tags")

(defvar *gff3-sequence-region-regex*
  (cl-ppcre:create-scanner
   "^##sequence-region\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)"))
(defvar *gff3-ontology-regex*
  (cl-ppcre:create-scanner "^##\\w+ontology\\s+(\\S+)"))

(defvar *parse-fn-symbols*
 '(parse-seqid parse-source parse-type
   parse-sequence-coord parse-sequence-coord
   parse-score parse-strand parse-phase
   parse-attributes)
  "GFF field parser function symcols, one for each of the 9 fields.")
(defvar *gff3-field-tags*
  '(:seqid :source :type :start :end
    :score :strand :phase :attributes)
  "GFF3 field keyword tags, one for each of the 9 fields.")

;; extra validation

;; CDS features must have a phase

;; ID attributes provided for features must be unique throughout the
;; gff3 file. These IDs are used to make part_of (Parent) associations
;; between features. Each line, if contains an ID, must be unique. For
;; multi-feature features, such as alignments, multiple lines can
;; share the same ID. If this is the case, the seqid, source, type
;; (method), strand, target name and all other attributes other than
;; target must be the same.

;; ignore ontology directives for now

;; unsolved issues:
;;
;; where in a vertex object do we put annotation?
;;

;; read a record; add it to the graph as a vertex; look for its
;; parent; if the parent is present, add a part_of relation; if the
;; parent is absent, defer relation by storing in a table with the
;; parent ID as key and child ID as value.


;; sequence-region is a bio-sequence (vertex)
;; gene is a bio-sequence (vertex) part_of sequence-region
;; mRNA is a bio-sequence (vertex) part_of gene

;; Make ontology-relationship edges, these have terms as predicates
;; predicate part_of
;; predicate derived_from


(defun wibble-gff3 (filename)
  (with-open-file (stream filename
                   :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((line-buffer (make-line-buffer stream))
          (graph (make-instance 'directed-acyclic-graph)))
      (unless (equalp (cons :gff-version 3)
                      (read-gff3 line-buffer))
        (error "does not appear to be GFF3 format~%"))
      (do ((gff (read-gff3 line-buffer) (read-gff3 line-buffer)))
          ((null gff) graph)
        (let ((record-type (assocdr :record-type gff)))
          (ecase record-type
            (:sequence-region
             (add-vertex (make-instance 'dna-sequence
                                        :identity
                                        (assocdr :seqid gff)
                                        :length
                                        (- (assocdr :end gff)
                                           (assocdr :start gff)))
                         graph))
            (:record
             (let ((attrs (assocdr :attributes gff)))
               (add-vertex (make-instance 'dna-sequence
                                          :identity
                                          (assocdr "ID" attrs :test #'equal))
                           graph)))))))))



(defun add-or-update-vertex (alist graph)
  (let* ((attrs (assocdr :attributes alist))
         (id (assocdr "ID" attrs :test #'equal))
         (vertex (lookup-vertex id graph)))
    (if vertex
        ;; (update-vertex vertex alist)
      (add-vertex (make-instance 'dna-sequence :identity id) graph))))

;; (defun update-vertex (vertex alist)
;;   (if (valid-vertex-update vertex alist)
;;       ))

(defun valid-vertex-update (vertex alist)
  "Checks that the :seqid :source :type :strand values in the update
ALIST agree with those already in VERTEX."
  (dolist (key '(:seqid :source :type :strand))
    (unless (equal (attribute-of vertex key)
                   (assocdr alist key))
      (error 'malformed-record-error :text
             (format nil "invalid ~a attribute (~a) in feature (~a); expected (~a) from previous record"
                     key (identity-of vertex) (assocdr alist key)
                     (attribute-of vertex key)))))
  t)


(defmethod read-gff3 ((line-buffer byte-line-buffer))
  (let ((line (find-line line-buffer #'content-bytes-p)))
    (when line
      (let ((str (make-sb-string line)))
        (cond ((gff3-directive-p str)
               (process-gff3-directive str))
              ((gff3-comment-p str)
               (process-gff3-comment str))
              (t
               (process-gff3-record str)))))))

(defun process-gff3-comment (str)
  "Returns an alist with key :comment and value being the entire
comment line STR, including the leading '#' character."
  (pairlis '(:record-type :comment) (list :comment str)))

(defun process-gff3-directive (str)
  "Returns an alist of parsed and tagged GFF directive data parsed
from STR."
  (cond ((string= "##gff-version" str :end2 13)
         (acons :record-type :gff-version
                (parse-gff-version str)))
        ((string= "##sequence-region" str :end2 17)
         (acons :record-type :sequence-region
                (parse-sequence-region str)))
        ((string= "##feature-ontology" str :end2 18)
         (acons :record-type :feature-ontology
                (parse-ontology-uri str)))
        ((string= "##attribute-ontology" str :end2 20)
         (acons :record-type :attribute-ontology
                (parse-ontology-uri str)))
        ((string= "##source-ontology" str :end2 17)
         (acons :record-type :source-ontology
                (parse-ontology-uri str)))
        ((string= "###" str :end2 3)
         (acons :record-type :end-forward-refs nil))
        ((string= "##FASTA" str :end2 7)
         (acons :record-type :fasta nil))
        (t
         (error 'malformed-record-error :text
                (format nil "unknown directive (~a)" str)))))

(defun process-gff3-record (str)
   "Returns an alist containing the record data. The alist keys are
:seqid :source :type :start :end :score :strand :phase and
:attributes. The value of :attributes is itself an alist, keyed by the
GFF attribute tag strings."
   (let ((record (parse-gff3-record str)))
     (unless (<= (assocdr :start record) (assocdr :end record))
       (error 'malformed-record-error :text
              (format nil "invalid feature coordinates (~a ~a); start must be <= end"
                      (assocdr :start record) (assocdr :end record))))
     record))

(defun parse-gff-version (str)
  "Returns an alist with key :version and value integer version
number parsed from STR. A value of 3 is expected."
  (handler-case
      (acons :version (parse-integer str :start 13) nil)
    (parse-error (condition)
      (error 'malformed-record-error :text (format nil "~a" condition)))))

(defun parse-sequence-region (str)
  "Returns an alist containing a seqid string, an integer sequence
start coordinate and an integer sequence end coordinate parsed from
STR."
  (multiple-value-bind (seqid start end)
      (cl-ppcre:register-groups-bind (x y z)
          (*gff3-sequence-region-regex* str)
        (values (parse-seqid x)
                (parse-sequence-coord y)
                (parse-sequence-coord z)))
    (unless (stringp seqid)
      (error 'malformed-record-error :text
             (format nil "invalid seqid (~a) in sequence-region directive (~a)"
                     seqid str)))
    (unless (integerp start)
      (error 'malformed-record-error :text
             (format nil "invalid start (~a) in sequence-region directive (~a)"
                     start str)))
    (unless (integerp end)
      (error 'malformed-record-error :text
             (format nil "invalid end (~a) in sequence-region directive (~a)"
                     end str)))
    (unless (<= start end)
      (error 'malformed-record-error :text
             (format nil "invalid sequence-region coordinates (~a ~a); start must be <= end"
                     start end)))
    (pairlis '(:seqid :start :end) (list seqid start end))))

(defun parse-gff3-record (str)
  "Returns alist of GFF record data parsed from STR."
  (multiple-value-bind (field-starts field-ends)
      (vector-split-indices #\Tab str)
    (unless (= 9 (length field-starts))
      (error 'malformed-record-error :text
             (format nil "invalid GFF line having ~a fields instead of 9 (~a)"
                     (length field-starts) str)))
    (let ((fields (mapcar #'(lambda (fn x y)
                              (funcall (symbol-function fn) str x y))
                          *parse-fn-symbols* field-starts field-ends)))
      (pairlis *gff3-field-tags* fields))))

(defun parse-ontology-uri (str)
  "Returns an alist with key :uri and value of an URI object parsed
from STR."
  (cl-ppcre:register-groups-bind (uri)
      (*gff3-ontology-regex* str)
    (acons :uri (puri:parse-uri uri) nil)))

(defun parse-seqid (str &optional (start 0) end)
  "Returns a seqid string extracted from line STR between START and
END."
  (let ((end (or end (length str))))
    (unless (loop for i from start below end
               always (or (gff3-valid-id-char-p (char str i))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "invalid seqid (~a)" (subseq str start end))))
    (subseq str start end)))

(defun parse-source (str &optional (start 0) end)
  "Returns a source string extracted from line STR between START and
END."
  (let ((end (or end (length str))))
    (unless (loop for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "invalid source (~a)" (subseq str start end))))
    (subseq str start end)))

(defun parse-type (str &optional (start 0) end)
  "Returns a type string extracted from line STR between START and
END."
  (let ((end (or end (length str))))
    (unless (loop for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "invalid type (~a)" (subseq str start end))))
    (subseq str start end)))

(defun parse-sequence-coord (str &optional (start 0) end)
  "Returns an integer sequence coordinate extracted from line STR
between START and END."
  (let* ((end (or end (length str)))
         (coord (handler-case
                    (parse-integer str :start start :end end)
                  (parse-error (condition)
                    (error 'malformed-record-error :text
                           (format nil "~a" condition))))))
    (unless (and (integerp coord)
                 (plusp coord))
      (error 'malformed-record-error :text
             (format nil "invalid sequence coordinate (~a); a
positive integer is required" coord)))
    coord))

(defun parse-score (str &optional (start 0) end)
  "Returns a float score from line STR between START and END. Returns
a float, or NIL if STR contains only the GFF empty column code between
START and END."
  (if (string= +gff3-dot-column+ str :start2 start :end2 end)
      nil
    (handler-case
        (parse-float str :start start :end end)
      (condition (condition)
        (when (subtypep (type-of condition) 'error)
          (error 'malformed-record-error :text
                 (format nil "~a" condition)))))))

(defun parse-strand (str &optional (start 0) end)
  "Returns a nucleic acid sequence strand object (canonical
SEQUENCE-STRAND instance) corresponding to line STR between START and
END."
  (cond ((string= +gff3-dot-column+ str :start2 start :end2 end)
         *without-strand*)
        ((string= +gff3-forward-strand+ str :start2 start :end2 end)
         *forward-strand*)
        ((string= +gff3-reverse-strand+ str :start2 start :end2 end)
         *reverse-strand*)
        ((string= +gff3-unknown-strand+ str :start2 start :end2 end)
         *unknown-strand*)
        (t
         (error 'malformed-record-error :text
                (format nil "invalid sequence strand (~a)"
                        (subseq str start end))))))

(defun parse-phase (str &optional (start 0) end)
  "Returns a coding phase integer from line STR between START and END,
or NIL if STR contains only the GFF empty column code between START
and END."
 (if (string= +gff3-dot-column+ str :start2 start :end2 end)
     nil
  (let ((phase (handler-case
                   (parse-integer str :start start :end end)
                 (parse-error (condition)
                   (error 'malformed-record-error :text
                          (format nil "~a" condition))))))
    (unless (and (integerp phase)
                 (<= 0 phase 2))
      (error 'malformed-record-error :text
             (format nil "invalid phase (~a); a positive integer between 0 and 2 (inclusive) is required"
                     phase)))
    phase)))

(defun parse-attributes (str &optional (start 0) end)
  "Returns an alist of GFF attributes from line STR between START and
END."
  (declare (optimize (speed 3) (debug 0)))
  (declare (type simple-string str))
  (let ((end (or end (length str))))
    (declare (type array-index start end))
    (unless (loop for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "invalid type (~a)" (subseq str start end))))
    (multiple-value-bind (attr-starts attr-ends)
        (vector-split-indices #\; str :start start :end end)
      (mapcar #'(lambda (x y)
                  (excise-attribute str x y))
              attr-starts attr-ends))))

(defun excise-attribute (str start end)
  "Returns a list containing a single GFF attribute extracted from
line STR between START and END. The first element of the list is the
key string and the rest of the list contains the value strings."
  (declare (optimize (speed 3) (debug 0)))
  (declare (type simple-string str)
           (type array-index start end))
  (let ((sep-index (position #\= str :start start :end end)))
    (cons (subseq str start sep-index)
          (vector-split #\, str :start (1+ sep-index) :end end))))

(defun content-bytes-p (bytes)
  "Returns T if BYTES does not consist entirely of whitespace
character codes, or NIL otherwise."
  (not (whitespace-bytes-p bytes)))

(defun gff3-directive-p (str)
  "Returns T if STR is a GFF3 directive line, or NIL otherwise."
  (and (>= (length str) 2)
       (char= #\# (char str 0) (char str 1))))

(defun gff3-comment-p (str)
  "Returns T if STR is a GFF3 comment line, or NIL otherwise."
  (and (> (length str) 0)
       (char= #\# (char str 0))
       (not (gff3-directive-p str))))

(defun gff3-valid-id-char-p (char)
  "Returns T if CHAR is a valid character for a GFF3 seqid, excluding
URL escapes."
  (loop for c in *gff3-valid-id-chars*
     thereis (char= char c)))

(defun url-escape-p (str index)
  "Returns T if INDEX into STR is the start of a valid URL encoded
character ('%' followed by two hexadecimal characters), or NIL
otherwise."
  (not (null (and (<= (+ index 3) (length str))
                  (char= #\% (char str index))
                  (digit-char-p (char str (+ index 1)) 16)
                  (digit-char-p (char str (+ index 2)) 16)))))
