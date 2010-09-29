;;;
;;; Copyright (C) 2010 Keith James. All rights reserved.
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

(defmethod hamming-search ((seq1 encoded-vector-sequence)
                           (seq2 encoded-vector-sequence)
                           &key (start1 0) end1 (start2 0) end2
                           (max-distance 1))
  (declare (optimize (speed 3) (safety 0)))
  (let ((vector1 (slot-value seq1 'vector))
        (vector2 (slot-value seq2 'vector)))
    (declare (type vector vector1 vector2)
             (type vector-index start1 start2)
             (type fixnum max-distance))
    (let* ((end1 (the vector-index (or end1 (length vector1))))
           (end2 (the vector-index (or end2 (length vector2))))
           (len1 (- end1 start1)))
      (loop
         with position = nil
         with distance = nil
         for i from start2 to (- end2 len1)
         for j of-type vector-index = (+ i len1)
         for d = (%hamming-distance vector1 vector2 :start1 start1 :end1 end1
                                    :start2 i :end2 j)
         for found = (<= d max-distance)
         do (when found
              (setf position i
                    distance d))
         until found
         finally (return (values position distance))))))

(defmethod hamming-distance ((seq1 encoded-vector-sequence)
                             (seq2 encoded-vector-sequence)
                             &key (start1 0) end1 (start2 0) end2)
  (let ((vector1 (slot-value seq1 'vector))
        (vector2 (slot-value seq2 'vector)))
    (let ((end1 (or end1 (length vector1)))
          (end2 (or end2 (length vector2))))
      (check-arguments (<= 0 start1 end1) (start1 end1)
                       "must satisfy (<= 0 start1 end1)")
      (check-arguments (<= 0 start2 end2) (start1 end2)
                       "must satisfy (<= 0 start2 end2)")
      ;; (check-arguments (= (- end1 start1)
      ;;                   (- end2 start2)) (start1 end1 start2 end2)
      ;;                   "ranges must be the same length, but were ~d and ~d"
      ;;                   (- end1 start1) (- end2 start2))
      (%hamming-distance vector1 vector2 :start1 start1 :end1 end1
                         :start2 start2 :end2 end2))))

(declaim (ftype (function (vector vector &key (:start1 fixnum) (:end1 t)
                                  (:start2 fixnum) (:end2 t)) fixnum)
                %hamming-distance))
(defun %hamming-distance (vector1 vector2
                          &key (start1 0) end1 (start2 0) end2)
  (macrolet ((hamming (vtype v1 s1 e1 v2 s2 e2)
               `(prog ()
                   (declare (optimize (speed 3) (safety 0)))
                   (declare (type ,vtype vector1 vector2))
                   (let ((,e1 (or ,e1 (length ,v1)))
                         (,e2 (or ,e2 (length ,v2))))
                     (return
                       (loop
                          for i of-type vector-index from ,s1 below ,e1
                          for j of-type vector-index from ,s2 below ,e2
                          when (/= (aref ,v1 i) (aref ,v2 j))
                          count i))))))
    ;; FIXME -- could add a clause for strings. Would need to
    ;; paramaterize the macrolet with a vector accessor and element
    ;; comparator
    (etypecase vector1
      ((encoded-residues 4)
       (hamming
        (encoded-residues 4) vector1 start1 end1 vector2 start2 end2))
      ((encoded-residues 7)
       (hamming
        (encoded-residues 7) vector1 start1 end1 vector2 start2 end2)))))
