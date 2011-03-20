;;;
;;; Copyright (c) 2009-2011 Keith James. All rights reserved.
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

;;; PowerLoom and Stella do not use ASDF. This system definition file
;;; may be copied to the root directory of a PowerLoom installation to
;;; allow loading of PowerLoom by ASDF.

;;; PowerLoom is a registered trademark of the University of Southern
;;; California.

(in-package :cl-user)

(defpackage :powerloom-system
  (:use :common-lisp :asdf))

(in-package :powerloom-system)

(defsystem powerloom
  :name "PowerLoom"
  :depends-on (:stella)
  :components ((:module :powerloom
                        :pathname ""
                        :components ((:file "load-powerloom")))))
