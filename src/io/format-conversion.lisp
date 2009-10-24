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

;; This is a naive converter that builds entire sequences in memory,
;; rather than streaming them
(defun convert-sequence-file (in-filespec in-format out-filespec out-format)
  "Converts the sequence data in the file identified
by IN-FILESPEC in format IN-FORMAT, to OUT-FORMAT, writing the data to
a new file identified by OUT-FILESPEC. Returns the number of records
converted."
  (with-ascii-li-stream (in in-filespec)
    (with-open-file (out out-filespec :direction :output
                         :element-type 'base-char
                         :external-format :ascii
                         :if-exists :supersede)
      (let ((gen (make-seq-input in in-format
                                 :parser (make-instance
                                          'raw-sequence-parser)))
            (con (make-seq-output out out-format)))
        (loop
           as alist = (next gen)
           while alist
           count alist
           do (funcall con alist))))))
