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

(defmacro in-syntax (readtable-expression)
  "Sets the readtable to the result of READTABLE-EXPRESSION for the
remainder of the file. Use within a file where the alternative
readtable is required. The original readtable will be restored once
the file is compiled or loaded.

This code by Kent Pitman, see X3J13 Cleanup Issue IN-SYNTAX:MINIMAL."
  `(eval-when (:compile-toplevel :load-toplevel :execute)
     (setq *readtable* ,readtable-expression)))

(defmacro define-pli-wrapper (fname (&rest args))
  "Defines a wrapper function which calls a PowerLoom function of the
same name. The wrapper function allows the module and environment
arguments required by many PLI functions to be made optional."
  `(defun ,fname (,@args &key (module *current-module*)
                  (env *current-environment*))
     (,(find-symbol (symbol-name fname) "PLI") ,@args module env)))

(defmacro with-environment ((env) &body body)
  "Evaluates BODY with the PowerLoom environment bound to ENV."
  `(let ((*current-environment* ,env))
    ,@body))

(defmacro with-module ((module) &body body)
  "Evaluates BODY with the PowerLoom module bound to MODULE."
  `(let ((*current-module* (if (symbolp ,module)
                               (get-module ,module)
                             ,module)))
    ,@body))

(defmacro in-module (module)
  "Sets *current-module* to PowerLoom MODULE. Analagous to CL:IN-PACKAGE."
  `(progn
     (setf *current-module* (if (symbolp ,module)
                                (get-module ,module)
                              module))
     (pli:change-module *current-module*)))

(defun dollar-reader (stream char)
  "This reader enables us to read PowerLoom concepts with dollar
syntax. Thus
;;; (get-concept 'concept-symbol)
becomes
;;; $concept-symbol
and
;;; (get-concept \"concept-name\")
becomes
;;; $\"concept-name\""
  (declare (ignore char))
  (let ((value (read stream t nil t)))
    (etypecase value
      (symbol `(get-concept ',value))
      (string `(s-get-concept ,value)))))

(defun setup-dollar-reader (&optional (readtable *readtable*))
  "Sets the {defun dollar-reader} to be used in READTABLE, dispatching
on the $ character."
  (set-macro-character #\$ #'dollar-reader nil readtable)
  readtable)

(defun ensure-string (name)
  (etypecase name
    (string name)
    (symbol (symbol-name name))))

(defun stella-symbol (name &optional (module *current-module*)
                      (env *current-environment*))
  "Returns a new Stella symbol given symbol or string NAME."
  (pli:create-symbol (ensure-string name) module env))

(defun to-stella (form &optional (module *current-module*)
                         (env *current-environment*))
  "Recursively converts FORM to a Stella equivalent. The result may be
used for PowerLoom queries."
  (cond ((symbolp form)
         (stella-symbol form module env))
        ((listp form)
         (mapcar (lambda (x)
                   (to-stella x module env)) form))
        ((integerp form)
         (stella::new-integer-wrapper form))
        ((floatp form)
         (stella::new-float-wrapper form))
        ((stringp form)
         (stella::new-string-wrapper form))
        (t
         form)))

(defun from-stella (object)
  "Returns the CL equivalent of Stella OBJECT."
  (cond ((pli:is-integer object)
         (pli:object-to-integer object))
        ((pli:is-float object)
         (pli:object-to-float object))
        (t
         (pli:object-to-string object))))

(defvar *powerloom-readtable*
  (setup-dollar-reader (copy-readtable nil))
  "The PowerLoom reader with dollar syntax for concepts.")
(defparameter *current-module* pli:null
  "The current PowerLoom module.")
(defparameter *current-environment* pli:null
  "The current PowerLoom environment.")

(define-pli-wrapper get-concept-instances (concept))
(define-pli-wrapper get-direct-concept-instances (concept))
(define-pli-wrapper get-direct-subrelations (relation))
(define-pli-wrapper get-direct-superrelations (relation))
(define-pli-wrapper get-direct-types (object))
(define-pli-wrapper is-subrelation (sub super))
(define-pli-wrapper is-a (object concept))

(defun load-file (filespec &optional (env *current-environment*))
  (etypecase filespec
    (string (pli:load filespec env))
    (pathname (pli:load (namestring filespec) env))
    (stream (pli:load-native-stream filespec env))))

(defun get-module (name &optional (env *current-environment*))
  (pli:get-module (ensure-string name) env))

(defun get-concept (name &key (module *current-module*)
                    (env *current-environment*))
  (pli:get-concept (ensure-string name) module env))

(defun evaluate (expression &key (module *current-module*)
                 (env *current-environment*))
  (pli:evaluate (to-stella expression module env) module env))

(defun ask (query &key (module *current-module*)
            (env *current-environment*))
  (pli:ask (to-stella query module env) module env))

(defun retrieve (query &key (module *current-module*)
                 (env *current-environment*))
  (pli:retrieve (to-stella query module env) module env))

(defun tuple (pl-iter &key (module *current-module*)
              (env *current-environment*))
  (loop
     for n from 0 below (pli:get-column-count pl-iter)
     collect (pli:get-nth-value pl-iter n module env)))

(defun tuple-size (pl-iter)
  (pli:get-column-count pl-iter))

(defun next-tuple (pl-iter)
  (when (pli:next? pl-iter)
    (tuple pl-iter)))

(defun collect-tuples (pl-iter)
  (loop
     while (pli:next? pl-iter)
     collect (tuple pl-iter)))

(defun nconc-tuples (pl-iter)
  (loop
     while (pli:next? pl-iter)
     nconc (tuple pl-iter)))

(defun subrelation-tree (relation &key (module *current-module*)
                         (env *current-environment*))
  (let ((subs (nconc-tuples (get-direct-subrelations
                             relation :module module :env env))))
    (if (null subs)
        relation
      (cons relation (mapcar (lambda (sub)
                               (subrelation-tree sub :module module :env env))
                             subs)))))
