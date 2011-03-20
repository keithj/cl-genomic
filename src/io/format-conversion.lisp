;;;
;;; Copyright (c) 2008-2011 Keith James. All rights reserved.
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

;; This is a naive converter that builds entire sequences in memory,
;; rather than streaming them
(defun convert-sequence-file (in-filespec in-format out-filespec out-format)
  "Converts the sequence data in the file identified
by IN-FILESPEC in format IN-FORMAT, to OUT-FORMAT, writing the data to
a new file identified by OUT-FILESPEC. Returns the number of records
converted."
  (with-seq-input (seqi in-filespec in-format
                        :parser (make-instance 'raw-sequence-parser))
    (with-seq-output (seqo out-filespec out-format)
      (loop
         for seq = (next seqi)
         while seq
         count seq
         do (consume seqo seq)))))
