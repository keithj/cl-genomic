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

(defmethod convert-sequence-file (in-filespec in-format
                                  out-filespec out-format)
  (error 'invalid-argument-error
         :params '(in-format out-format)
         :args (list in-format out-format)
         :text "automatic conversion not available"))

(defmethod convert-sequence-file (in-filespec (in-format (eql :fastq))
                                  out-filespec (out-format (eql :fasta)))
  (with-ascii-li-stream (in in-filespec)
    (with-open-file (out out-filespec :direction :output
                         :element-type 'base-char
                         :external-format :ascii
                         :if-exists :supersede)
      (let ((gen (make-seq-input in in-format
                                 :parser (make-instance
                                          'raw-sequence-parser))))
        (loop
           as raw = (next gen)
           while raw
           count raw
           do (write-raw-fasta raw out))))))

(defmethod convert-sequence-file (in-filespec (in-format (eql :fasta))
                                  out-filespec (out-format (eql :fasta)))
  (with-ascii-li-stream (in in-filespec)
    (with-open-file (out out-filespec :direction :output
                         :element-type 'base-char
                         :external-format :ascii
                         :if-exists :supersede)
      (let ((gen (make-seq-input in in-format
                                 :parser (make-instance
                                          'raw-sequence-parser))))
        (loop
           as raw = (next gen)
           while raw
           count raw
           do (write-raw-fasta raw out))))))
