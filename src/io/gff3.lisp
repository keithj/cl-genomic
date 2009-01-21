;;;
;;; Copyright (C) 2008-2009 Keith James. All rights reserved.
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

;; (defmethod read-gff3 ((stream character-line-input-stream)
;;                       (parser gff3-parser))
;;   (loop
;;      for line = (find-line stream #'content-string-p)
;;      while (not (eql :eof line))
;;      do (cond ((gff3-directive-p line)
;;                (gff3-directive parser (parse-gff3-directive line)))
;;               ((gff3-comment-p line)
;;                (gff3-comment parser (parse-gff3-comment line)))
;;               (t
;;                (gff3-record parser (parse-gff3-record line))))))




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
        (error 'malformed-record-error
               :record new
               :text (format nil "invalid record: conflicts in fields ~a"
                             conflicts)))
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

(defun parse-gff3-comment (str)
  "Returns an alist with key :comment and value being the entire
comment line STR, including the leading '#' character."
  (pairlis '(:record-type :comment) (list :comment str)))

(defun parse-gff3-directive (str)
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
         (error 'malformed-record-error
                :record str
                :text "unknown directive"))))

(defun parse-gff3-record (str)
   "Returns an alist containing the record data. The alist keys are
:seqid :source :type :start :end :score :strand :phase and
:attributes. The value of :attributes is itself an alist, keyed by the
GFF attribute tag strings."
   (multiple-value-bind (field-starts field-ends)
       (vector-split-indices #\Tab str)
     (unless (= (length *gff3-field-tags*) (length field-starts))
       (error 'malformed-record-error
             :record str
             :text (format nil "invalid record having ~a fields instead of ~a"
                           (length field-starts)
                           (length *gff3-field-tags*))))
     (let ((record (pairlis
                    *gff3-field-tags*
                    (mapcar (lambda (parse-fn x y)
                              (funcall parse-fn str x y))
                            (field-parse-fns) field-starts field-ends))))
       (unless (<= (assocdr :start record) (assocdr :end record))
         (error 'malformed-record-error
                :record record
                :text "invalid feature coordinates having start > end"))
       record)))

(defun parse-gff-version (str)
  "Returns an alist with key :version and value integer version
number parsed from STR. A value of 3 is expected."
  (handler-case
      (acons :version (parse-integer str :start 13) nil)
    (parse-error (condition)
      (error 'malformed-record-error
             :record str
             :text (format nil "~a" condition)))))

(defun parse-sequence-region (str)
  "Returns an alist containing a seqid string, an integer sequence
start coordinate and an integer sequence end coordinate parsed from
STR."
  (cl-ppcre:register-groups-bind (x y z)
      (*gff3-sequence-region-regex* str)
    (let ((seqid (parse-seqid x))
          (start (parse-gff3-sequence-coord y))
          (end (parse-gff3-sequence-coord z)))
      (let ((record (pairlis '(:seqid :start :end) (list seqid start end))))
        (unless (<= start end)
          (error 'malformed-record-error
                 :record record
                 :text "invalid sequence-region coordinates having start > end"))
        record))))

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
      (error 'malformed-field-error
             :record str
             :field (subseq str start end)
             :text "invalid seqid"))
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
      (error 'malformed-field-error
             :record str
             :field (subseq str start end)
             :text "invalid source"))
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
      (error 'malformed-field-error
             :record str
             :field (subseq str start end)
             :text "invalid type"))
    (subseq str start end)))

(defun parse-gff3-sequence-coord (str &optional (start 0) end)
  "Returns an integer sequence coordinate extracted from line STR
between START and END."
  (let* ((end (or end (length str)))
         (coord (handler-case
                    (parse-integer str :start start :end end)
                  (parse-error (condition)
                    (error 'malformed-record-error
                           :text (format nil "~a" condition))))))
    (unless (and (integerp coord)
                 (plusp coord))
      (error 'malformed-field-error
             :record str
             :field coord
             :text "invalid sequence coordinate: a positive integer is required"))
    coord))

(defun parse-gff3-score (str &optional (start 0) end)
  "Returns a float score from line STR between START and END. Returns
a float, or NIL if STR contains only the GFF empty column code between
START and END."
  (if (string= +gff3-dot-column+ str :start2 start :end2 end)
      nil
    (handler-case
        (parse-float str :start start :end end)
      (error ()
        (error 'malformed-field-error
               :record str
               :field (subseq str start end)
               :text "invalid score")))))

(defun parse-gff3-strand (str &optional (start 0) end)
  "Returns a nucleic acid sequence strand object (canonical
SEQUENCE-STRAND instance) corresponding to line STR between START and
END."
  (if (string= +gff3-dot-column+ str :start2 start :end2 end)
      nil ; FIXME -- was *without-strand*, should it be nil?
    (decode-strand str :strict t)))

(defun parse-gff3-phase (str &optional (start 0) end)
  "Returns a coding phase integer from line STR between START and END,
or NIL if STR contains only the GFF empty column code between START
and END."
 (if (string= +gff3-dot-column+ str :start2 start :end2 end)
     nil
  (let ((phase (handler-case
                   (parse-integer str :start start :end end)
                 (parse-error ()
                   (error 'malformed-field-error
                          :record str
                          :field (subseq str start end)
                          :text "invalid phase")))))
    (unless (and (integerp phase)
                 (<= 0 phase 2))
      (error 'malformed-field-error
             :record str
             :field phase
             :text "invalid phase: a positive integer between 0 and 2 (inclusive) is required"))
    phase)))

(defun parse-gff3-attributes (str &optional (start 0) end)
  "Returns an alist of GFF attributes from line STR between START and
END."
  (let ((end (or end (length str))))
    (unless (loop
               for i from start below end
               always (or (not (control-char-p (char str i)))
                          (when (char= #\% (char str i))
                            (url-escape-p str i))))
      (error 'malformed-record-error
             :record str
             :text "unescaped control characters in attributes"))
    (multiple-value-bind (attr-starts attr-ends)
        (vector-split-indices #\; str :start start :end end)
      (mapcar (lambda (x y)
                (excise-attribute str x y)) attr-starts attr-ends))))

(defun excise-attribute (str start end)
  "Returns a list containing a single GFF attribute extracted from
line STR between START and END. The first element of the list is the
key string and the rest of the list contains the value strings."
  (let ((sep-index (position #\= str :start start :end end)))
    (cons (subseq str start sep-index)
          (vector-split #\, str :start (1+ sep-index) :end end))))

(defun gff-version-p (record version)
  "Returns T if RECORD is a GFF version directive (:record-type
:gff-version) indicating VERSION. If RECORD is not a GFF version
directive an error is thrown."
  (unless (eql :gff-version (assocdr :record-type record))
    (error 'malformed-record-error
           :record record
           :text "invalid GFF version directive"))
  (eql version (assocdr :version record)))

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

(defun field-parse-fns ()
  "Returns a list of GFF field parser functions, one for each of the 9
fields."
  (list #'parse-seqid #'parse-gff3-source #'parse-gff3-type
        #'parse-gff3-sequence-coord #'parse-gff3-sequence-coord
        #'parse-gff3-score #'parse-gff3-strand #'parse-gff3-phase
        #'parse-gff3-attributes))
