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

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defvar *codons*
    '((#\t #\t #\t) (#\t #\t #\c) (#\t #\t #\a) (#\t #\t #\g) (#\t #\c #\t)
      (#\t #\c #\c) (#\t #\c #\a) (#\t #\c #\g) (#\t #\a #\t) (#\t #\a #\c)
      (#\t #\a #\a) (#\t #\a #\g) (#\t #\g #\t) (#\t #\g #\c) (#\t #\g #\a)
      (#\t #\g #\g) (#\c #\t #\t) (#\c #\t #\c) (#\c #\t #\a) (#\c #\t #\g)
      (#\c #\c #\t) (#\c #\c #\c) (#\c #\c #\a) (#\c #\c #\g) (#\c #\a #\t)
      (#\c #\a #\c) (#\c #\a #\a) (#\c #\a #\g) (#\c #\g #\t) (#\c #\g #\c)
      (#\c #\g #\a) (#\c #\g #\g) (#\a #\t #\t) (#\a #\t #\c) (#\a #\t #\a)
      (#\a #\t #\g) (#\a #\c #\t) (#\a #\c #\c) (#\a #\c #\a) (#\a #\c #\g)
      (#\a #\a #\t) (#\a #\a #\c) (#\a #\a #\a) (#\a #\a #\g) (#\a #\g #\t)
      (#\a #\g #\c) (#\a #\g #\a) (#\a #\g #\g) (#\g #\t #\t) (#\g #\t #\c)
      (#\g #\t #\a) (#\g #\t #\g) (#\g #\c #\t) (#\g #\c #\c) (#\g #\c #\a)
      (#\g #\c #\g) (#\g #\a #\t) (#\g #\a #\c) (#\g #\a #\a) (#\g #\a #\g)
      (#\g #\g #\t) (#\g #\g #\c) (#\g #\g #\a) (#\g #\g #\g))
    "The codons, ordered for making genetic codes."))

(defvar *genetic-codes* (make-hash-table)
  "The standard genetic codes.")

(defmacro define-genetic-code (symbol &key name full-name id
                               amino-acids starts
                               (aa-encoder #'encode-aa-7bit)
                               (na-encoder #'encode-dna-4bit))
  "Defines a new genetic code and registers it.

Arguments:

- symbol (symbol): A symbol used to name a dynamic variable to which
  the new genetic code object is bound.

Key:

- name (symbol): The genetic code name, e.g. :standard
- full-name (string): The genetic code full name e.g. \"Standard\"
- id (integer): The genetic code number.
- amino-acids (string): A 64 element string of amino acid tokens.
- starts (string): A 64 element string of amino acid tokens indicating
  potential start codons.
- aa-encoder (function): A function used to encode the amino acid
  tokens.
- na-encoder (function): A function used to encode the codon tokens."
  (let ((t1 (make-hash-table :test #'equal))
        (t2 (make-hash-table :test #'equal))
        (t3 (make-hash-table :test #'equal))
        (asp-asn (mapcar #'encode-aa-7bit '(#\D #\N)))  ; D and N -> B
        (glu-gln (mapcar #'encode-aa-7bit '(#\E #\Q)))) ; E and Q -> Z
    (mapc (lambda (codon aa start)
            (let ((encoded-codon (mapcar na-encoder codon))
                  (encoded-aa (funcall aa-encoder aa)))
              (setf (gethash encoded-codon t1) encoded-aa)
              (when (char= #\M start)
                (setf (gethash encoded-codon t2) t))
              (when (char= #\* aa)
                (setf (gethash encoded-codon t3) t))))
          *codons* (coerce amino-acids 'list) (coerce starts 'list))
    (with-gensyms (trans-table start-table term-table)
      `(progn
        (defvar ,symbol (make-instance 'genetic-code
                                       :name ,name :full-name ,full-name
                                       :identity ,id))
        (register-genetic-code ,symbol)
        (let ((,trans-table ,t1)
              (,start-table ,t2)
              (,term-table ,t3))
          (defmethod translate-codon (codon (code (eql ,symbol))
                                      &key initiator)
            (flet ((trn (c)
                     (gethash c ,trans-table)))
              (let ((encoded-aa
                     (delete-duplicates
                      (mapcar #'trn (enum-encoded-codon codon)))))
                (cond ((= 1 (length encoded-aa))
                       (cond ((and initiator (gethash codon ,start-table))
                              (encode-aa-7bit #\M))
                             (initiator
                              (error 'initiator-codon-error
                                     :codon (mapcar #'decode-rna-4bit codon)
                                     :genetic-code code))
                             (t
                              (first encoded-aa))))
                      ((= 2 (length encoded-aa))
                       (cond ((null (set-difference ',glu-gln encoded-aa))
                              (encode-aa-7bit #\B)) ; D and N -> B
                             ((null (set-difference ',asp-asn encoded-aa))
                              (encode-aa-7bit #\Z)) ; E and Q -> Z
                             (t
                              (encode-aa-7bit #\X))))
                      (t
                       (encode-aa-7bit #\X))))))
          (defmethod start-codon-p (codon (code (eql ,symbol)))
            (gethash codon ,start-table))
          (defmethod term-codon-p (codon (code (eql ,symbol)))
            (gethash codon ,term-table)))))))

(defun register-genetic-code (genetic-code)
  "Registers a global standard GENETIC-CODE."
  (setf (gethash (name-of genetic-code) *genetic-codes*) genetic-code))

(defun find-genetic-code (name)
  "Returns a standard GENETIC-CODE designated by its name symbol NAME,
such as :standard, :vert-mito or :bacterial."
  (multiple-value-bind (genetic-code presentp)
      (gethash name *genetic-codes*)
    (unless presentp
      (error 'invalid-argument-error
             :params 'name
             :args name
             :text "no such genetic code"))
    genetic-code))

(defun registered-genetic-codes ()
  "Returns a list of all registered genetic codes."
  (loop
     for code being the hash-values of *genetic-codes*
     collect code into codes
     finally (return (sort codes #'< :key #'identity-of))))

(defclass genetic-code (identity-mixin)
  ((name :initform nil
         :initarg :name
         :reader name-of
         :documentation "The genetic code symbolic name.")
   (full-name :initform nil
              :initarg :full-name
              :reader full-name-of
              :documentation "The genetic code full name."))
  (:documentation "A genetic code translating codons to amino acids."))

(defmethod print-object ((code genetic-code) stream)
  (format stream "#<GENETIC-CODE ~a ~a>" (slot-value code 'identity)
          (slot-value code 'name)))

(define-genetic-code *standard*
    :name :standard :id 1
    :full-name "Standard"
    :amino-acids
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "---M---------------M---------------M----------------------------")

(define-genetic-code *vert-mito*
    :name :vert-mito :id 2
    :full-name "Vertebrate Mitochondrial"
    :amino-acids
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
    :starts
    "--------------------------------MMMM---------------M------------")

(define-genetic-code *yeast-mito*
    :name :yeast-mito :id 3
    :full-name "Yeast Mitochondrial"
    :amino-acids
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "----------------------------------MM----------------------------")

(define-genetic-code *mold-mito*
    :name :mold-mito :id 4
    :full-name "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"
    :amino-acids
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "--MM---------------M------------MMMM---------------M------------")

(define-genetic-code *invert-mito*
    :name :invert-mito :id 5
    :full-name "Invertebrate Mitochondrial"
    :amino-acids
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
    :starts
    "---M----------------------------MMMM---------------M------------")

(define-genetic-code *ciliate-nuclear*
    :name :ciliate-nuclear :id 6
    :full-name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"
    :amino-acids
    "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M----------------------------")

;;; No genetic codes 7 and 8

(define-genetic-code *echino-mito*
    :name :echino-mito :id 9
    :full-name "Echinoderm Mitochondrial; Flatworm Mitochondrial"
    :amino-acids
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M---------------M------------")

(define-genetic-code *euplotid-nuc*
    :name :euplotid-nuc :id 10
    :full-name "Euplotid Nuclear"
    :amino-acids
    "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M----------------------------")

(define-genetic-code *bacterial*
    :name :bacterial :id 11
    :full-name "Bacterial and Plant Plastid"
    :amino-acids
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "---M---------------M------------MMMM---------------M------------")

(define-genetic-code *alt-yeast-nuc*
    :name :alt-yeast-nuc :id 12
    :full-name "Alternative Yeast Nuclear"
    :amino-acids
    "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "-------------------M---------------M----------------------------")

(define-genetic-code *ascidian-nuc*
    :name :alt-ascidian-nuc :id 13
    :full-name "Ascidian Mitochondrial"
    :amino-acids
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
    :starts
    "---M------------------------------MM---------------M------------")

(define-genetic-code *alt-flatworm-nuc*
    :name :alt-flatworm-nuc :id 14
    :full-name "Alternative Flatworm Mitochondrial"
    :amino-acids
    "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M----------------------------")

(define-genetic-code *blepharisma-macronuc*
    :name :blepharisma-macronuc :id 15
    :full-name "Blepharisma Macronuclear"
    :amino-acids
    "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M----------------------------")

(define-genetic-code *chlorophycean-mito*
    :name :chlorophycean-mito :id 16
    :full-name "Chlorophycean Mitochondrial"
    :amino-acids
    "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M----------------------------")

;;; No genetic codes 17, 18, 19 and 20

(define-genetic-code *trematode-mito*
    :name :trematode-mito :id 21
    :full-name "Trematode Mitochondrial"
    :amino-acids
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M---------------M------------")

(define-genetic-code *scenedesmus-mito*
    :name :scenedesmus-mito :id 22
    :full-name "Scenedesmus obliquus Mitochondrial"
    :amino-acids
    "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "-----------------------------------M----------------------------")

(define-genetic-code *thraustochytrium-mito*
    :name :thraustochytrium-mito :id 23
    :full-name "Thraustochytrium Mitochondrial"
    :amino-acids
    "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    :starts
    "--------------------------------M--M---------------M------------")
