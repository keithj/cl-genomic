;;;
;;; Copyright (C) 2008-2009 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
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

(defvar *default-base-concept* nil)

(defvar *obo-escape-chars*
  '(#\n #\W #\t #\: #\, #\" #\\ #\( #\) #\[ #\] #\{ #\})
  "Characters that may escaped in OBO tags and values.")

(defvar *obo-header-mandatory-tags*
  '("format-version"))

(defvar *obo-header-optional-tags*
  '("auto-generated-by" "date" "saved-by" "version" "data-version"
    "default-relationship-id-prefix" "id-mapping" "idspace" "remark"
    "subsetdef" "synonymtypedef"))

(defvar *obo-builtin-identifers*
  '("OBO:TYPE" "OBO:TERM" "OBO:TERM_OR_TYPE" "OBO:INSTANCE"))

(defvar *obo-builtin-primitives*
  '("xsd:boolean" "xsd:date" "xsd:decimal" "xsd:integer"
    "xsd:negativeInteger" "xsd:nonNegativeInteger"
    "xsd:nonPositiveInteger" "xsd:positiveInteger"
    "xsd:string" "xsd:simpleType"))

(defvar *obo-term-mandatory-tags*
  '("id" "name"))

(defvar *obo-term-optional-tags*
  '("alt_id" "builtin" "comment" "consider" "def" "disjoint_from"
    "intersection_of" "is_a" "is_obsolete" "relationship" "replaced_by"
    "subset" "synonym" "union_of" "use_term" "xref" "is_anonymous"))

(defvar *obo-term-deprecated-tags*
  '("exact_synonym" "narrow_synonym" "xref_analog" "xref_unk"
    "broad_synonym"))

(defvar *obo-typedef-mandatory-tags*
  '("id" "name"))

(defvar *obo-typedef-optional-tags*
  (sort (union
         (set-difference *obo-term-optional-tags*
                         '("disjoint_from" "intersection_of" "union_of"))
         '("domain" "range" "inverse_of" "is_anti_symmetric" "is_cyclic"
           "is_metadata_tag" "is_reflexive" "is_symmetric" "is_transitive"
           "transitive_over")) #'string<))

(defvar *obo-instance-mandatory-tags*
  '("id" "instance_of" "name"))

(defvar *obo-expected-format-version* "1.2")

(defun read-obo-stream (stream parser)
  (loop
     as line = (read-wrapped-line stream)
     while (not (eql :eof line))
     do (cond ((comment-line-p line)
               nil)
              ((whitespace-string-p line)
               (end-section parser))
              ((term-stanza-p line)
               (begin-term parser))
              ((typedef-stanza-p line)
               (begin-typedef parser))
              ((instance-stanza-p line)
               (begin-instance parser))
              (t
               (multiple-value-bind (tag value)
                   (read-tag-value line)
                 (tag-value parser tag value))))
     finally (end-section parser))
  parser)

;;; Default methods which do nothing except collect tag-values and
;;; check that the mandatory ones are present
(defmethod begin-term ((parser obo-parser))
  nil)

(defmethod begin-typedef ((parser obo-parser))
  nil)

(defmethod begin-instance ((parser obo-parser))
  nil)

(defmethod end-section ((parser obo-parser))
  nil)

(defmethod begin-term :before ((parser obo-parser))
  (with-accessors ((state state-of) (tag-values tag-values-of))
      parser
    (setf state 'term
          tag-values ())))

(defmethod begin-typedef :before ((parser obo-parser))
  (with-accessors ((state state-of) (tag-values tag-values-of))
      parser
    (setf state 'typedef
          tag-values ())))

(defmethod begin-instance :before ((parser obo-parser))
  (with-accessors ((state state-of) (tag-values tag-values-of))
      parser
    (setf state 'instance
          tag-values ())))

(defmethod end-section :before ((parser obo-parser))
  (with-accessors ((state state-of) (tag-values tag-values-of))
      parser
    (ecase state
      (header (check-mandatory-tags tag-values
                                    *obo-header-mandatory-tags*))
      (term (check-mandatory-tags tag-values
                                  *obo-term-mandatory-tags*)
            (check-tag-counts tag-values))
      (typedef (check-mandatory-tags tag-values
                                     *obo-typedef-mandatory-tags*)
               (check-tag-counts tag-values))
      (instance (check-mandatory-tags tag-values
                                      *obo-instance-mandatory-tags*)
                (check-tag-counts tag-values))
      ((nil) nil)))) ; case for multiple or trailing empty lines

(defmethod tag-value ((parser obo-parser) tag value)
  (with-accessors ((tag-values tag-values-of))
      parser
    (setf tag-values (acons tag value tag-values))))

(defun comment-line-p (str)
  "Returns T if STR is a comment line, or NIL otherwise."
  (starts-with-char-p (string-left-trim '(#\Space) str) #\!))

(defun term-stanza-p (str)
  "Returns T if STR is the header for a term stanza, or NIL
otherwise."
  (starts-with-string-p str "[Term]"))

(defun typedef-stanza-p (str)
  "Returns T if STR is the header for a typedef stanza, or NIL
otherwise."
  (starts-with-string-p str "[Typedef]"))

(defun instance-stanza-p (str)
  "Returns T if STR is the header for an instance stanza, or NIL
otherwise."
  (starts-with-string-p str "[Instance]"))

(defun read-tag-value (str)
  "Parses STR into tag and value strings, returning them as two
values. The tags and values retain escape characters, trailing
modifiers, dbxref lists and comments."
  (let ((i (unescaped-position #\: str))
        (trim '(#\Space)))
    (values (string-trim trim (subseq str 0 i))
            (string-trim trim (subseq str (1+ i))))))

(defun read-value (tag-value)
  (let ((tag (car tag-value))
        (value (cdr tag-value)))
    (cons tag (cond ((string= "name" tag)
                     (read-name value))
                    ((string= "def" tag)
                     (read-def value))
                    ((string= "is_a" tag)
                     (read-is-a value))
                    ((string= "relationship" tag)
                     (read-relationship value))
                    ((string= "intersection_of" tag)
                     (read-intersection value))
                    (t
                     value)))))


;;; Functions for specific tags
(defun read-def (str)
  "Returns the quoted text portion of STR. Currently ignores any
dbxref list."
  (read-quoted-text str))

(defun read-is-a (str)
  (expand-escape-chars (remove-trailing-modifiers (remove-comments str))))

(defun read-name (str)
  (expand-escape-chars (remove-comments str)))

(defun read-relationship (str)
  (split-value str))

(defun read-intersection (str)
  (split-value str))

(defun split-value (str)
  (let ((parts (string-split (remove-comments str) #\Space
                             :remove-empty-substrings t)))
    (if (endp (rest parts))
        (string-trim '(#\Space) (first parts))
      parts)))

;;; Functions for all tags
(defun read-quoted-text (str)
  (let ((start (unescaped-position #\" str)))
    (cond ((null start)
           nil)
          ((= (length str) (1+ start))
           (error 'malformed-field-error
                  :field str
                  :text "unbalanced quotes"))
          (t
           (let ((end (unescaped-position #\" str :start (1+ start))))
             (if end
                 (subseq str (1+ start) end)
               (error 'malformed-field-error
                  :field str
                  :text "unbalanced quotes")))))))

(defun remove-comments (str)
  (let ((excl-index (position #\! str)))
    (if excl-index
        (subseq str 0 excl-index)
      str)))

(defun remove-trailing-modifiers (str)
  (let ((brace-index (unescaped-position #\{ str :from-end t)))
    (string-trim '(#\Space) (if brace-index
                                (subseq str 0 brace-index)
                              str))))

(defun remove-dbxref-list (str)
  (let ((bracket-index (unescaped-position #\[ str  :from-end t)))
    (string-trim '(#\Space) (if bracket-index
                                (subseq str 0 bracket-index)
                              str))))

(defun unescaped-position (char str &key (start 0) end from-end)
  "Returns the position of the first unescaped CHAR in STR, between
START (which defaults to 0) and END, or NIL otherwise."
  (let* ((end (or end (length str)))
         (count (- end start)))
    (loop
       for i in (if from-end
                    (iota count (1- end) -1)
                  (iota count start))
       do (when (and (char= char (char str i))
                     (not (escape-char-p str i)))
            (return i)))))

(defun wrapped-line-p (str)
  "Returns T if STR ends with an escaped literal newline, indicating a
wrapped line, or NIL otherwise."
  (ends-with-char-p str #\\))

(defun read-wrapped-line (stream)
  "Reads a line, which may be wrapped, from STREAM."
  (loop
     for line = (stream-read-line stream)
     until (eql :eof line)
     collect line into lines
     while (wrapped-line-p line)
     finally (return (if (eql :eof line)
                         :eof
                       (join-wrapped-lines lines)))))

(defun join-wrapped-lines (lines)
  "Returns a string created by joining the list of strings with
escaped newlines LINES."
  (flet ((remove-escape (line)
           (subseq line 0  (1- (length line)))))
    (if (null lines)
        nil
      (apply #'concatenate 'string
             (append (mapcar #'remove-escape (butlast lines))
                     (last lines))))))

(defun escape-char-p (str index)
  "Returns T if there is an OBO escape character at INDEX in STR, or
NIL otherwise."
  (or (and (char= #\\ (char str index))
           (< (1+ index) (length str))
           (find (char str (1+ index)) *obo-escape-chars*))
      (and (find (char str index) *obo-escape-chars*)
           (< 0 index (length str))
           (char= #\\ (char str (1- index))))))

(defun expand-escape-chars (str)
  "Returns a new string created by unescaping all the OBO escape
characters in STR."
  (do ((expanded ())
       (len (length str))
       (i 0 (1+ i)))
      ((= i len) (make-array (length expanded)
                             :element-type 'character
                             :initial-contents (nreverse expanded)))
    (cond ((escape-char-p str i)
           (incf i)
           (push (expand-escape-char (char str i)) expanded))
          (t
           (push (char str i) expanded)))))

(defun expand-escape-char (char)
  "Returns the character represented by the OBO escape character
CHAR."
  (case char
    (#\n #\Newline)
    (#\W #\Space)
    (#\t #\Tab)
    (otherwise char)))

(defun check-mandatory-tags (tag-values mandatory-tags)
  (loop
     for tag in mandatory-tags
     do (unless (assocdr tag tag-values :test #'string=)
          (error 'malformed-record-error
                 :record tag-values
                 :text (format nil "missing mandatory tag ~a" tag))))
  t)

(defun check-id-value (tag-values)
  (let ((id (assocdr "id" tag-values :test #'string=)))
    (loop
       for tag in (append *obo-builtin-identifers*
                          *obo-builtin-primitives*)
       do (when (string= id tag)
            (error 'malformed-record-error
                   :record tag-values
                   :text (format nil "~a is not a valid id tag" id)))))
  t)

(defun check-tag-counts (tag-values)
  (loop
     for (tag . value) in tag-values  ; how to declare that VALUE is ignored?
     count (string= "id" tag) into id-count
     count (string= "def" tag) into def-count
     count (string= "name" tag) into name-count
     count (string= "comment" tag) into comment-count
     finally (cond ((> id-count 1)
                    (error 'malformed-record-error
                           :record tag-values
                           :text ">1 id tag present"))
                   ((> def-count 1)
                    (error 'malformed-record-error
                           :record tag-values
                           :text ">1 def tag present"))
                   ((> name-count 1)
                    (error 'malformed-record-error
                           :record tag-values
                           :text ">1 name tag present"))
                   ((> comment-count 1)
                    (error 'malformed-record-error
                           :record tag-values
                           :text ">1 comment tag present"))))
  t)
