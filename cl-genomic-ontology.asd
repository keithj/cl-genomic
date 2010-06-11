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

(asdf:load-system :deoxybyte-systems)

(in-package :deoxybyte-systems)

(defsystem cl-genomic-ontology
  :depends-on (:cl-genomic)
  :in-order-to ((test-op (load-op :cl-genomic-ontology
                                  :cl-genomic-ontology-test)))
  :components ((:module :powerloom
                        :serial t
                        :pathname "src/ontology/"
                        :components
                        ((:file "powerloom-loader")
                         (:file "package")
                         (:file "powerloom")
                         (:file "sofa-knowledgebase")))
               (:lift-test-config :lift-tests
                                  :pathname "cl-genomic-ontology-test"
                                  :target-system :cl-genomic-ontology)))
