;;;
;;; Copyright (C) 2008 Keith. All rights reserved.
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

(defstruct (interval)
  "A biological sequence interval bounded by a pair of inter-residue
coordinates."
  (lower 0 :type fixnum)
  (upper 0 :type fixnum))

(defstruct (stranded-interval (:include interval)
                              (:conc-name s-interval-))
  "A biological sequence interval specific to a sequence strand."
  (strand *unknown-strand* :type sequence-strand))

