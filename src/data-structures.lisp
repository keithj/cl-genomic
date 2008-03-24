;;;
;;; Copyright (C) 2008 Keith James. All rights reserved.
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

(defun make-seq-datum (identity alphabet ambiguity &key token-seq
                       length description)
  "Returns a biological sequence data structure, given a sequence
IDENTITY, a vector of residue tokens TOKEN-SEQ and a DESCRIPTION
string. The intention is that constructing one of these is cheaper
than making a CLOS instance from the same data. This is useful for
cases such as round-trip IO."
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

(defun make-quality-datum (identity alphabet ambiguity &key token-seq
                           length quality)
  "Returns a biological sequence with quality data structure, given a
sequence IDENTITY, a vector of residue tokens TOKEN-SEQ and a QUALITY
vector."
  (acons :quality quality
         (make-seq-datum identity alphabet ambiguity
                         :token-seq token-seq :length length)))

(defun seq-datum-identity (datum)
  "Returns the biological sequence identity from DATUM."
  (assocdr :identity datum))

(defun seq-datum-alphabet (datum)
  "Returns the biological sequence alphabet from DATUM."
  (assocdr :alphabet datum))

(defun seq-datum-ambiguity (datum)
  "Returns the biological sequence ambiguity from DATUM."
  (assocdr :ambiguity datum))

(defun seq-datum-token-seq (datum)
  "Returns the biological sequence token-seq from DATUM."
  (assocdr :token-seq datum))

(defun seq-datum-length (datum)
  "Returns the biological sequence length from DATUM."
  (assocdr :length datum))

(defun seq-datum-description (datum)
  "Returns the biological sequence description from DATUM."
  (assocdr :description datum))

(defun seq-datum-quality (datum)
  "Returns the biological sequence quality from DATUM."
  (assocdr :quality datum))

(defun make-seq-from-datum (datum)
  "A callback which constructs a CLOS bio-sequence object from a
sequence datum."
  (make-seq :alphabet (seq-datum-alphabet datum)
            :ambiguity (seq-datum-ambiguity datum)
            :identity (seq-datum-identity datum)
            :token-seq (seq-datum-token-seq datum)
            :length (seq-datum-length datum)))
