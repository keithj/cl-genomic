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

(defmethod read-bio-sequence :before (stream format &key alphabet ambiguity
                                             virtualp)
  (unless alphabet
    (error "An alphabet must be supplied."))
  (unless (member ambiguity '(:iupac :auto nil))
    (error "Invalid ambiguity ~a: expected one of ~a."
           ambiguity '(:iupac :auto nil)))
  (unless (member virtualp '(t nil))
    (error "Invalid virtualp ~a: expected one of ~a."
           virtualp '(t nil))))


(defun make-seq-alist (identity alphabet ambiguity &key token-seq
                       length description)
  "Returns an alist, given a sequence IDENTITY, a vector of residue
tokens TOKEN-SEQ and a DESCRIPTION string."
  (pairlis '(:identity :alphabet :ambiguity :token-seq :length :description)
           (list identity alphabet
                 (cond ((and (eql :auto ambiguity)
                             (not (simplep token-seq alphabet)))
                        :iupac)
                       ((eql :auto ambiguity)
                        nil)
                       (t
                        ambiguity))
                 token-seq length description)))

(defun make-quality-alist (identity alphabet ambiguity &key token-seq
                           length quality)
  "Returns an alist, given a sequence IDENTITY, a vector of residue
tokens TOKEN-SEQ and a QUALITY vector."
  (acons :quality quality
         (make-seq-alist identity alphabet ambiguity
                         :token-seq token-seq :length length)))

(defun make-seq-from-alist (alist)
  "A callback which constructs a CLOS bio-sequence object from
sequence data that has been parsed into an ALIST."
  (make-seq  :alphabet (assocdr :alphabet alist)
             :ambiguity (assocdr :ambiguity alist)
             :identity (assocdr :identity alist)
             :token-seq (assocdr :token-seq alist)
             :length (assocdr :length alist)))

(defun make-chunk-pname (file-pname chunk-number)
  "Returns a new pathname for a file chunk based on a file pathname
FILE-PNAME and an integer CHUNK-NUMBER."
  (make-pathname
   :directory (pathname-directory file-pname)
   :name (concatenate 'string (pathname-name file-pname) "."
                      (princ-to-string chunk-number))
   :type (pathname-type file-pname)))
