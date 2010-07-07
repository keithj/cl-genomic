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
  (:setup (unless (modulep :sequence-ontology)
            (load-ontology (merge-pathnames "ontology/so_2_4_3.plm"))
            (load-ontology (merge-pathnames "ontology/so_addenda.plm")))))

(in-syntax *powerloom-readtable*)

 (addtest (cl-genomic-ontology-tests) modulep/1
   (ensure (modulep "SEQUENCE-ONTOLOGY"))
   (ensure (modulep :sequence-ontology))
   (ensure (not (modulep :foo))))

 (addtest (cl-genomic-ontology-tests) get-module/1
   (ensure (get-module "SEQUENCE-ONTOLOGY"))
   (ensure (get-module :sequence-ontology)))

 (addtest (cl-genomic-ontology-tests) with-module/1
   (ensure (equal pli:null *current-module*))
   (with-module (:sequence-ontology)
     (ensure (equal *current-module* (get-module :sequence-ontology))))
   (ensure (equal pli:null *current-module*)))

 (addtest (cl-genomic-ontology-tests) ampersand-reader/1
  (with-module (:sequence-ontology)
    (ensure (is-logic-object @SO:0000110))))

(addtest (cl-genomic-ontology-tests) term-name/1
  (with-module (:sequence-ontology)
    (ensure (equal "sequence_feature" (object-to-string (term-name @SO:0000110))))))

(addtest (cl-genomic-ontology-tests) term-doc/1
  (with-module (:sequence-ontology)
    (ensure (is-string (term-doc @SO:0000110)))))

(addtest (cl-genomic-ontology-tests) find-term/1
  (with-module (:sequence-ontology)
    (ensure (equal @SO:0000110 (find-term "sequence_feature")))))

(addtest (cl-genomic-ontology-tests) find-doc/1
  (with-module (:sequence-ontology)
    (ensure (is-string (find-doc "sequence_feature")))
    (ensure (equal (term-doc @SO:0000110) (find-doc "sequence_feature")))))

(addtest (cl-genomic-ontology-tests) term-parents/1
  (with-module (:sequence-ontology)
    ;; gene (SO:0000704) has parents gene_group (via part_of) and
    ;; biological_region, region and sequence_feature (via is_a)
    (ensure (dxu:linear-set-equal
             (list @SO:0005855 @SO:0001411 @SO:0000001 @SO:0000110)
             (term-parents @SO:0000704) :test #'equal :key #'get-name))))

(addtest (cl-genomic-ontology-tests) term-parent-p/1
  (with-module (:sequence-ontology)
    ;; gene (SO:0000704) has parents gene_group (via part_of) and
    ;; biological_region, region and sequence_feature (via is_a)
    (dolist (term (list @SO:0005855 @SO:0001411 @SO:0000001 @SO:0000110))
      (ensure (term-parent-p term @SO:0000704)))))

