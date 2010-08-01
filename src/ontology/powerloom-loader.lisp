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

(in-package :cl-user)

(defvar *powerloom-home* nil
  "Path to the PowerLoom installation.")
(defvar *powerloom-loader* "load-powerloom.lisp"
  "The name of the loader file supplied with PowerLoom.")

(defun load-powerloom ()
  (flet ((load-pl (pl-home)
           (handler-bind (#+sbcl(sb-ext:compiler-note #'muffle-warning)
                                (style-warning #'muffle-warning)
                                (warning #'muffle-warning))
             ;; Swallow PowerLoom's chatter
             (let ((*standard-output* (make-broadcast-stream)))
               (load (merge-pathnames *powerloom-loader*
                                      (fad:pathname-as-directory pl-home)))))))
    (cond (*powerloom-home*
           (load-pl *powerloom-home*))
          ((dxi:environment-variable "POWERLOOM_HOME")
           (load-pl (dxi:environment-variable "POWERLOOM_HOME")))
          (t
           (error "The *powerloom-home* was not set. See README.txt for help.")))))

;; These are setq'd in PowerLoom's load-stella file, but are not
;; defined on SBCL. Defining them here avoids compiler warnings.
#+:sbcl(defparameter *stella-compiler-optimization* nil)
#+:sbcl(defparameter *gc-verbose* nil)

(load-powerloom)
(funcall (intern "INITIALIZE" 'pli))
