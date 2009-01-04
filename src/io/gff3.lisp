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
  "GFF field parser function symbols, one for each of the 9 fields.")
(defvar *gff3-field-tags*
  '(:seqid :source :type :start :end
    :score :strand :phase :attributes)
  "GFF3 field keyword tags, one for each of the 9 fields.")
(defvar *gff3-invariant-field-tags*
  '(:seqid :source :type :score :strand))

(defvar *gff3-reserved-attribute-tags*
  '("ID" "Name" "Alias" "Parent" "Target" "Gap" "Derives_from"
    "Note" "Dbxref" "Ontology_term" "Index")
  "Reserved GFF3 attribute tags")
(defvar *gff3-invariant-attribute-tags*
  '("ID" "Name" "Alias" "Parent" "Gap" "Derives_from"
    "Note" "Dbxref" "Ontology_term" "Index"))

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

(defparameter *expected-gff-version* 3)

;; (defun wibble-gff3 (filename)
;;   (with-open-file (fs filename
;;                    :direction :input
;;                    :element-type '(unsigned-byte 8))
;;     (let ((stream (make-line-input-stream fs))
;;           (graph (make-instance 'directed-acyclic-graph))
;;           (records (make-hash-table :test #'equal)))
;;       (unless (gff-version-p (read-gff3 stream) *expected-gff-version*)
;;         (error "Does not appear to be GFF3 format."))
;;       (do ((record (read-gff3 stream) (read-gff3 stream)))
;;           ((null record) (make-sequence-features records graph))
;;         (let ((record-type (assocdr :record-type record)))
;;           (ecase record-type
;;             (:sequence-region
;;              (make-sequence-region record graph))
;;             (:record
;;              (open-or-update-record record records))
;;             (:end-forward-refs
;;              (make-sequence-features records graph))))))))

;; (defun make-sequence-region (record graph)
;;   "Returns a new sequence-region vertex created from the data in
;; RECORD. A a side-effect, adds the sequence-region to GRAPH."
;;   (when (lookup-vertex (assocdr :seqid record) graph)
;;     (error 'malformed-record-error :text
;;            (format nil (msg "Invalid sequence-region ~a:"
;;                             "a region with this seqid has already"
;;                             "been found."))))
;;   (add-vertex (make-instance 'dna-sequence
;;                              :identity (assocdr :seqid record)
;;                              :length (- (assocdr :end record)
;;                                         (assocdr :start record)))
;;   )

(defun open-or-update-record (record records)
  "Adds the data in RECORD to the table of currently open RECORDS."
  (let ((identity (assocdr "ID" (assocdr :attributes record)
                           :test #'equal)))
    (setf (gethash identity records)
          (merge-records (gethash identity records) record))))

(defun merge-records (new current)
  (if (null current)
      (copy-alist new)
    (let ((conflicts (loop
                        for field in *gff3-invariant-field-tags*
                        when (not (eql (assocdr field new)
                                       (assocdr field current)))
                        collect field)))
      (when conflicts
        (error 'malformed-record-error :text
               "Invalid record: conflicts in fields ~a." conflicts))
      (loop
         for field in (set-difference *gff3-field-tags*
                                      *gff3-invariant-field-tags*)
         do (assocpush+ field current (assocdr field new)))
      ;; FIXME -- deal with attributes here
      )))

;; (defun valid-vertex-update (vertex alist)
;;   "Checks that the :seqid :source :type :strand values in the update
;; ALIST agree with those already in VERTEX."
;;   (dolist (key '(:seqid :source :type :strand))
;;     (unless (equal (attribute-of vertex key)
;;                    (assocdr alist key))
;;       (error 'malformed-record-error :text
;;              (format nil (msg "Invalid ~a attribute ~a in feature ~a:"
;;                               "expected ~a from previous record.")
;;                      key (identity-of vertex) (assocdr alist key)
;;                      (attribute-of vertex key)))))
;;   t)


(defun gff-version-p (record version)
  "Returns T if RECORD is a GFF version directive (:record-type
:gff-version) indicating VERSION. If RECORD is not a GFF version
directive an error is thrown."
  (unless (eql :gff-version (assocdr :record-type record))
    (error "Record ~a is not a GFF version directive." record))
  (eql version (assocdr :version record)))





(defmethod read-gff3 ((stream binary-line-input-stream))
  (let ((line (find-line stream #'content-bytes-p)))
    (if (vectorp line)
        (let ((str (make-sb-string line)))
          (cond ((gff3-directive-p str)
                 (process-gff3-directive str))
                ((gff3-comment-p str)
                 (process-gff3-comment str))
                (t
                 (process-gff3-record str))))
      nil)))

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
                (format nil "Unknown directive ~a." str)))))

(defun process-gff3-record (str)
   "Returns an alist containing the record data. The alist keys are
:seqid :source :type :start :end :score :strand :phase and
:attributes. The value of :attributes is itself an alist, keyed by the
GFF attribute tag strings."
   (let ((record (parse-gff3-record str)))
     (unless (<= (assocdr :start record) (assocdr :end record))
       (error 'malformed-record-error :text
              (format nil (msg "Invalid feature coordinates (~a ~a):"
                               "start must be <= end.")
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
                (parse-gff3-sequence-coord y)
                (parse-gff3-sequence-coord z)))
    (unless (stringp seqid)
      (error 'malformed-record-error :text
             (format nil "Invalid seqid ~a in sequence-region directive ~a."
                     seqid str)))
    (unless (integerp start)
      (error 'malformed-record-error :text
             (format nil "Invalid start ~a in sequence-region directive ~a."
                     start str)))
    (unless (integerp end)
      (error 'malformed-record-error :text
             (format nil "Invalid end ~a in sequence-region directive ~a."
                     end str)))
    (unless (<= start end)
      (error 'malformed-record-error :text
             (format nil (msg "Invalid sequence-region coordinates (~a ~a)"
                              "start must be <= end.")
                     start end)))
    (pairlis '(:seqid :start :end) (list seqid start end))))

(defun parse-gff3-record (str)
  "Returns alist of GFF record data parsed from STR."
  (multiple-value-bind (field-starts field-ends)
      (vector-split-indices #\Tab str)
    (unless (= 9 (length field-starts))
      (error 'malformed-record-error :text
             (format nil (msg "Invalid GFF line having ~a fields"
                              "instead of 9: (~a).")
                     (length field-starts) str)))
    (let ((fields (mapcar (lambda (fn x y)
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
    (unless (loop
               for i from start below end
               always (or (gff3-valid-id-char-p (char str i))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "Invalid seqid ~a." (subseq str start end))))
    (subseq str start end)))

(defun parse-gff3-source (str &optional (start 0) end)
  "Returns a source string extracted from line STR between START and
END."
  (let ((end (or end (length str))))
    (unless (loop
               for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "Invalid source ~a." (subseq str start end))))
    (subseq str start end)))

(defun parse-gff3-type (str &optional (start 0) end)
  "Returns a type string extracted from line STR between START and
END."
  (let ((end (or end (length str))))
    (unless (loop
               for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "Invalid type ~a." (subseq str start end))))
    (subseq str start end)))

(defun parse-gff3-sequence-coord (str &optional (start 0) end)
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
             (format nil (msg "Invalid sequence coordinate ~a:"
                              "a positive integer is required.")
                     coord)))
    coord))

(defun parse-gff3-score (str &optional (start 0) end)
  "Returns a float score from line STR between START and END. Returns
a float, or NIL if STR contains only the GFF empty column code between
START and END."
  (if (string= +gff3-dot-column+ str :start2 start :end2 end)
      nil
    (handler-case
        (parse-float str :start start :end end)
      (condition (condition)
        (when (subtypep (type-of condition) 'error) ; what other conditions?
          (error 'malformed-record-error :text
                 (format nil "Invalid score ~a." condition)))))))

(defun parse-gff3-strand (str &optional (start 0) end)
  "Returns a nucleic acid sequence strand object (canonical
SEQUENCE-STRAND instance) corresponding to line STR between START and
END."
  (cond ((string= +gff3-dot-column+ str :start2 start :end2 end)
         nil) ; FIXME -- was *without-strand*, should it be nil?
        ((string= +gff3-forward-strand+ str :start2 start :end2 end)
         *forward-strand*)
        ((string= +gff3-reverse-strand+ str :start2 start :end2 end)
         *reverse-strand*)
        ((string= +gff3-unknown-strand+ str :start2 start :end2 end)
         *unknown-strand*)
        (t
         (error 'malformed-record-error :text
                (format nil "Invalid sequence strand ~a."
                        (subseq str start end))))))

(defun parse-gff3-phase (str &optional (start 0) end)
  "Returns a coding phase integer from line STR between START and END,
or NIL if STR contains only the GFF empty column code between START
and END."
 (if (string= +gff3-dot-column+ str :start2 start :end2 end)
     nil
  (let ((phase (handler-case
                   (parse-integer str :start start :end end)
                 (parse-error (condition)
                   (error 'malformed-record-error :text
                          (format nil "Invalid phase ~a." condition))))))
    (unless (and (integerp phase)
                 (<= 0 phase 2))
      (error 'malformed-record-error :text
             (format nil (msg "Invalid phase ~a: a positive integer"
                              "between 0 and 2 (inclusive) is required.")
                     phase)))
    phase)))

(defun parse-gff3-attributes (str &optional (start 0) end)
  "Returns an alist of GFF attributes from line STR between START and
END."
  (declare (optimize (speed 3) (debug 0)))
  (declare (type simple-string str))
  (let ((end (or end (length str))))
    (declare (type array-index start end))
    (unless (loop
               for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error :text
             (format nil "Invalid type ~a." (subseq str start end))))
    (multiple-value-bind (attr-starts attr-ends)
        (vector-split-indices #\; str :start start :end end)
      (mapcar (lambda (x y)
                (excise-attribute str x y)) attr-starts attr-ends))))

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
  (not (null (and (< (+ index 2) (length str))
                  (char= #\% (char str index))
                  (digit-char-p (char str (+ index 1)) 16)
                  (digit-char-p (char str (+ index 2)) 16)))))
