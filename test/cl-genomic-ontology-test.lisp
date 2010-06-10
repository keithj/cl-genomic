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

(in-package :cl-genomic-ontology-test)

;; The base test of all cl-genomic-ontology tests
(deftestsuite cl-genomic-ontology-tests ()
  ()
  (:setup (unless (modulep :sofa)
            (load-ontology (merge-pathnames "ontology/sofa_2_4_2.plm"))
            (load-ontology (merge-pathnames "ontology/sofa_addenda.plm")))))

(addtest (cl-genomic-ontology-tests) modulep/1
  (ensure (modulep "SOFA"))
  (ensure (modulep :sofa))
  (ensure (not (modulep :foo))))

(addtest (cl-genomic-ontology-tests) get-module/1
  (ensure (get-module "SOFA"))
  (ensure (get-module :sofa)))
