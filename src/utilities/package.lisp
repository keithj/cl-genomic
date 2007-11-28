
(defpackage bigvector
  (:use #:common-lisp)
  (:documentation "Specialised vectors of unlimited size.")
  (:export
   ;; Functions
   ;; Classes
   #:bigvector
   ;; Generic functions
   #:length-of
   #:element-type-of
   #:bvref
   #:bigvector-subseq
   #:bigvector-concatenate
   #:bigvector-copy
   ;; Conditions
   #:illegal-subscript
   #:incompatible-start-and-end
   #:subscript-out-of-bounds))