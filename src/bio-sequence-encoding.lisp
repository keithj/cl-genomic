;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
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
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See tnhe
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :bio-sequence)

(defun complement-dna (base)
  "Returns the complentary DNA base to that represented by the
character BASE. Folds to lower-case."
  (ecase base
    ((#\t #\T) #\a)
    ((#\c #\C) #\g)
    ((#\a #\A) #\t)
    ((#\g #\G) #\c)
    ((#\r #\R) #\y)
    ((#\y #\Y) #\r)
    ((#\k #\K) #\m)
    ((#\m #\M) #\k)
    ((#\s #\S) #\s)
    ((#\w #\W) #\w)
    ((#\b #\B) #\v)
    ((#\d #\D) #\h)
    ((#\h #\H) #\d)
    ((#\v #\V) #\b)
    ((#\n #\N) #\n)
    (#\-       #\-)))

(defun complement-rna (base)
  "Returns the complentary RNA base to that represented by the
character BASE. Folds to lower-case."
  (ecase base
    ((#\u #\U) #\a)
    ((#\c #\C) #\g)
    ((#\a #\A) #\u)
    ((#\g #\G) #\c)
    ((#\r #\R) #\y)
    ((#\y #\Y) #\r)
    ((#\k #\K) #\m)
    ((#\m #\M) #\k)
    ((#\s #\S) #\s)
    ((#\w #\W) #\w)
    ((#\b #\B) #\v)
    ((#\d #\D) #\h)
    ((#\h #\H) #\d)
    ((#\v #\V) #\b)
    ((#\n #\N) #\n)
    (#\-       #\-)))

(declaim (inline complement-dna-8bit))
(defun complement-dna-8bit (byte)
  "Returns the byte representing the complementary DNA base to that of
BYTE. Folds to lower case."
  (ecase byte
    ((116 84)  97)
    ((99  67) 103)
    ((97  65) 116)
    ((103 71)  99)
    ((114 82) 121)
    ((121 89) 114)
    ((107 75) 109)
    ((109 77) 107)
    ((115 83) 115)
    ((119 87) 119)
    ((98  66) 118)
    ((100 68) 104)
    ((104 72) 100)
    ((118 86)  98)
    ((110 78) 110)
    (45        45)))

(declaim (inline complement-dna-4bit))
(defun complement-dna-4bit (encoded-base)
  "Returns the complentary encoded DNA base to that represented by
ENCODED-BASE."
  (ecase encoded-base
    (#b0001 #b1000)
    (#b0010 #b0100)
    (#b0100 #b0010)
    (#b1000 #b0001)
    (#b0101 #b1010)                     ; purine
    (#b1010 #b0101)                     ; pyrimidine
    (#b1100 #b0011)                     ; keto
    (#b0011 #b1100)                     ; amino
    (#b0110 #b0110)                     ; strong
    (#b1001 #b1001)                     ; weak
    (#b1110 #b0111)                     ; !adenosine
    (#b1101 #b1011)                     ; !cytosine
    (#b1011 #b1101)                     ; !guanine
    (#b0111 #b1110)                     ; !thymine
    (#b1111 #b1111)                     ; any
    (#b0000 #b0000)))                   ; gap

(declaim (inline complement-rna-4bit))
(defun complement-rna-4bit (encoded-base)
  "Returns the complentary encoded RNA base to that represented by
ENCODED-BASE."
   (ecase encoded-base
    (#b0001 #b1000)
    (#b0010 #b0100)
    (#b0100 #b0010)
    (#b1000 #b0001)
    (#b0101 #b1010)                     ; purine
    (#b1010 #b0101)                     ; pyrimidine
    (#b1100 #b0011)                     ; keto
    (#b0011 #b1100)                     ; amino
    (#b0110 #b0110)                     ; strong
    (#b1001 #b1001)                     ; weak
    (#b1110 #b0111)                     ; !adenosine
    (#b1101 #b1011)                     ; !cytosine
    (#b1011 #b1101)                     ; !guanine
    (#b0111 #b1110)                     ; !thymine
    (#b1111 #b1111)                     ; any
    (#b0000 #b0000)))                   ; gap

(declaim (inline encode-dna-4bit))
(defun encode-dna-4bit (base)
  "Encodes DNA standard-char BASE as a 4-bit byte, representing A as
0001, C as 0010, G as 0100 and T as 1000. The first base is in the
most significant 4-bit byte and the last base is in the least
significant 4-bit byte. Ambiguous bases are represented by bitwise AND
combinations of these."
  (ecase base
    ((#\a #\A) #b0001)
    ((#\c #\C) #b0010)
    ((#\g #\G) #b0100)
    ((#\t #\T) #b1000)
    ((#\r #\R) #b0101)                  ; purine
    ((#\y #\Y) #b1010)                  ; pyrimidine
    ((#\k #\K) #b1100)                  ; keto
    ((#\m #\M) #b0011)                  ; amino
    ((#\s #\S) #b0110)                  ; strong
    ((#\w #\W) #b1001)                  ; weak
    ((#\b #\B) #b1110)                  ; !adenosine
    ((#\d #\D) #b1101)                  ; !cytosine
    ((#\h #\H) #b1011)                  ; !guanine
    ((#\v #\V) #b0111)                  ; !thymine
    ((#\n #\N) #b1111)                  ; any
    (#\-       #b0000)))                ; gap

(declaim (inline decode-dna-4bit))
(defun decode-dna-4bit (encoded-base)
  "Decodes 4-bit byte ENCODED-BASE, returning a lower case character."
  (ecase encoded-base
    (#b0001 #\a)
    (#b0010 #\c)
    (#b0100 #\g)
    (#b1000 #\t)
    (#b0101 #\r)
    (#b1010 #\y)
    (#b1100 #\k)
    (#b0011 #\m)
    (#b0110 #\s)
    (#b1001 #\w)
    (#b1110 #\b)
    (#b1101 #\d)
    (#b1011 #\h)
    (#b0111 #\v)
    (#b1111 #\n)
    (#b0000 #\-)))

(declaim (inline encode-dna-comp-4bit))
(defun encode-dna-comp-4bit (base)
  "Encodes the complement of DNA standard-char BASE as a 4-bit byte,
representing T as 0001, C as 0010, A as 0100 and G as 1000. The first
base is in the most significant 4-bit byte and the last base is in the
least significant 4-bit byte. Ambiguous bases are represented by
bitwise AND combinations of these."
  (encode-dna-4bit (complement-dna base)))

(declaim (inline decode-dna-comp-4bit))
(defun decode-dna-comp-4bit (encoded-base)
   "Decodes 4-bit byte ENCODED-BASE, returning its complement as a
lower case character."
  (complement-dna (decode-dna-4bit encoded-base)))

(declaim (inline encode-rna-4bit))
(defun encode-rna-4bit (base)
  "Encodes RNA standard-char BASE as a 4-bit byte, representing U as
0001, C as 0010, A as 0100 and G as 1000.  The first base is in the
most significant 4-bit byte and the last base is in the least
significant 4-bit byte. Ambiguous bases are represented by bitwise AND
combinations of these."
  (ecase base
    ((#\a #\A) #b0001)
    ((#\c #\C) #b0010)
    ((#\g #\G) #b0100)
    ((#\u #\U) #b1000)
    ((#\r #\R) #b0101)                  ; purine
    ((#\y #\Y) #b1010)                  ; pyrimidine
    ((#\k #\K) #b1100)                  ; keto
    ((#\m #\M) #b0011)                  ; amino
    ((#\s #\S) #b0110)                  ; strong
    ((#\w #\W) #b1001)                  ; weak
    ((#\b #\B) #b1110)                  ; !adenosine
    ((#\d #\D) #b1101)                  ; !cytosine
    ((#\h #\H) #b1011)                  ; !guanine
    ((#\v #\V) #b0111)                  ; !thymine
    ((#\n #\N) #b1111)                  ; any
    (#\-       #b0000)))                ; gap

(declaim (inline decode-rna-4bit))
(defun decode-rna-4bit (encoded-base)
  "Decodes 4-bit byte ENCODED-BASE, returning a lower case character."
  (ecase encoded-base
    (#b0001 #\a)
    (#b0010 #\c)
    (#b0100 #\g)
    (#b1000 #\u)
    (#b0101 #\r)
    (#b1010 #\y)
    (#b1100 #\k)
    (#b0011 #\m)
    (#b0110 #\s)
    (#b1001 #\w)
    (#b1110 #\b)
    (#b1101 #\d)
    (#b1011 #\h)
    (#b0111 #\v)
    (#b1111 #\n)
    (#b0000 #\-)))

(declaim (inline encode-rna-comp-4bit))
(defun encode-rna-comp-4bit (base)
  "Encodes the complement of RNA standard-char BASE as a 4-bit byte,
representing U as 0001, C as 0010, A as 0100 and G as 1000.  The first
base is in the most significant 4-bit byte and the last base is in the
least significant 4-bit byte. Ambiguous bases are represented by
bitwise AND combinations of these."
  (encode-rna-4bit (complement-rna base)))

(declaim (inline decode-rna-comp-4bit))
(defun decode-rna-comp-4bit (encoded-base)
  "Decodes 4-bit byte ENCODED-BASE, returning its complement as a
lower case character."
  (complement-rna (decode-rna-4bit encoded-base)))

(declaim (inline encode-aa-7bit))
(defun encode-aa-7bit (aa)
  "Encodes AA standard-char as a 7-bit byte."
  (ccase aa
    ((#\a #\A) 1)
    ((#\b #\B) 2)
    ((#\c #\C) 3)
    ((#\d #\D) 4)
    ((#\e #\E) 5)
    ((#\f #\F) 6)
    ((#\g #\G) 7)
    ((#\h #\H) 8)
    ((#\i #\I) 9)
    ((#\j #\J) 10)
    ((#\k #\K) 11)
    ((#\l #\L) 12)
    ((#\m #\M) 13)
    ((#\n #\N) 14)
    ((#\o #\O) 15)
    ((#\p #\P) 16)
    ((#\q #\Q) 17)
    ((#\r #\R) 18)
    ((#\s #\S) 19)
    ((#\t #\T) 20)
    ((#\u #\U) 21)
    ((#\v #\V) 22)
    ((#\w #\W) 23)
    ((#\x #\X) 24)
    ((#\y #\Y) 25)
    ((#\z #\Z) 26)
    (#\*       27)
    (#\-       0)))

(declaim (inline decode-aa-7bit))
(defun decode-aa-7bit (encoded-aa)
  "Decodes 7-bit byte ENCODED-AA, returning an upper case character."
  (ecase encoded-aa
    (1  #\A)
    (2  #\B)
    (3  #\C)
    (4  #\D)
    (5  #\E)
    (6  #\F)
    (7  #\G)
    (8  #\H)
    (9  #\I)
    (10 #\J)
    (11 #\K)
    (12 #\L)
    (13 #\M)
    (14 #\N)
    (15 #\O)
    (16 #\P)
    (17 #\Q)
    (18 #\R)
    (19 #\S)
    (20 #\T)
    (21 #\U)
    (22 #\V)
    (23 #\W)
    (24 #\X)
    (25 #\Y)
    (26 #\Z)
    (27 #\*)
    (0  #\-)))

(defun symbolize-dna-base (base)
  "Returns a Lisp keyword symbol representing the base character
BASE."
  (let ((encoded-base (encode-dna-4bit base)))
    (ecase encoded-base
      (#b0001 :thymine)
      (#b0010 :cytosine)
      (#b0100 :adenine)
      (#b1000 :guanine)
      (#b0101 :purine)
      (#b1010 :pyrimidine)
      (#b1100 :keto)
      (#b0011 :amino)
      (#b0110 :strong)
      (#b1001 :weak)
      (#b1110 :!adenine)
      (#b1101 :!cytosine)
      (#b1011 :!guanine)
      (#b0111 :!thymine)
      (#b1111 :any)
      (#b0000 :gap))))

(defun encode-dna-symbol (dna-symbol)
  "Encodes the symbol DNA-SYMBOL as a 4-bit byte."
  (ecase dna-symbol
    (:thymine    #b0001)
    (:cytosine   #b0010)
    (:adenine    #b0100)
    (:guanine    #b1000)
    (:purine     #b0101)
    (:pyrimidine #b1010)
    (:keto       #b1100)
    (:amino      #b0011)
    (:strong     #b0110)
    (:weak       #b1001)
    (:!adenine   #b1110)
    (:!cytosine  #b1101)
    (:!guanine   #b1011)
    (:!thymine   #b0111)
    (:any        #b1111)
    (:gap        #b0000)))

(defun symbolize-rna-base (base)
  "Returns a Lisp keyword symbol representing the base character
BASE."
  (let ((encoded-base (encode-rna-4bit base)))
    (ecase encoded-base
      (#b0001 :uracil)
      (#b0010 :cytosine)
      (#b0100 :adenine)
      (#b1000 :guanine)
      (#b0101 :purine)
      (#b1010 :pyrimidine)
      (#b1100 :keto)
      (#b0011 :amino)
      (#b0110 :strong)
      (#b1001 :weak)
      (#b1110 :!adenine)
      (#b1101 :!cytosine)
      (#b1011 :!guanine)
      (#b0111 :!thymine)
      (#b1111 :any)
      (#b0000 :gap))))

(defun encode-rna-symbol (rna-symbol)
  "Encodes the symbol RNA-SYMBOL as a 4-bit byte."
  (ecase rna-symbol
    (:uracil     #b0001)
    (:cytosine   #b0010)
    (:adenine    #b0100)
    (:guanine    #b1000)
    (:purine     #b0101)
    (:pyrimidine #b1010)
    (:keto       #b1100)
    (:amino      #b0011)
    (:strong     #b0110)
    (:weak       #b1001)
    (:!adenine   #b1110)
    (:!cytosine  #b1101)
    (:!guanine   #b1011)
    (:!thymine   #b0111)
    (:any        #b1111)
    (:gap        #b0000)))

(defun symbolize-aa (aa)
  "Returns a Lisp keyword symbol representing the amino-acid character
AA."
  (let ((encoded-aa (encode-aa-7bit aa)))
    (ecase encoded-aa
      (1  :alanine)
      (2  :aspartatic-acid/asparagine)
      (3  :cysteine)
      (4  :aspartic-acid)
      (5  :glutamic-acid)
      (6  :phenylalanine)
      (7  :glycine)
      (8  :histidine)
      (9  :isoleucine)
      (10 :isoleucine/leucine)
      (11 :lysine)
      (12 :leucine)
      (13 :methionine)
      (14 :asparagine)
      (15 :pyrrolysine)
      (16 :proline)
      (17 :glutamine)
      (18 :arginine)
      (19 :serine)
      (20 :threonine)
      (21 :selenocysteine)
      (22 :valine)
      (23 :tryptophan)
      (24 :unknown)
      (25 :tyrosine)
      (26 :glutamic-acid/glutamine)
      (27 :terminator)
      (0  :gap))))

(defun encode-aa-symbol (aa-symbol)
  "Encodes the symbol AA-SYMBOL as a 7-bit byte."
  (ecase aa-symbol
    (:alanine 1)
    (:aspartatic-acid/asparagine 2)
    (:cysteine 3)
    (:aspartic-acid 4)
    (:glutamic-acid 5)
    (:phenylalanine 6)
    (:glycine 7)
    (:histidine 8)
    (:isoleucine 9)
    (:isoleucine/leucine 10)
    (:lysine 11)
    (:leucine 12)
    (:methionine 13)
    (:asparagine 14)
    (:pyrrolysine 15)
    (:proline 16)
    (:glutamine 17)
    (:arginine 18)
    (:serine 19)
    (:threonine 20)
    (:selenocysteine 21)
    (:valine 22)
    (:tryptophan 23)
    (:unknown 24)
    (:tyrosine 25)
    (:glutamic-acid/glutamine 26)
    (:terminator 27)
    (:gap 0)))

(defun enum-dna-base (base)
  "Returns a sorted list of the unambiguous DNA characters represented
by DNA character BASE."
  (sort (mapcar #'decode-dna-4bit
                (enum-encoded-base (encode-dna-4bit base)))
        #'char<=))

(defun enum-rna-base (base)
  "Returns a sorted list of the unambiguous RNA characters represented
by RNA character BASE."
  (sort (mapcar #'decode-rna-4bit
                (enum-encoded-base (encode-rna-4bit base)))
        #'char<=))

(defun enum-dna-codon (codon)
  "Returns a list of the unambiguous DNA character codons represented
by list of DNA characters CODON."
  (mapcar (lambda (ecodon)
            (mapcar #'decode-dna-4bit ecodon))
          (enum-encoded-codon
           (mapcar #'encode-dna-4bit codon))))

(defun enum-rna-codon (codon)
  "Returns a list of the unambiguous RNA character codons represented
by list of RNA characters CODON."
  (mapcar (lambda (ecodon)
            (mapcar #'decode-rna-4bit ecodon))
          (enum-encoded-codon
           (mapcar #'encode-rna-4bit codon))))

(defun enum-aa (aa)
  "Returns a list of the unambiguous AA characters represented by
amino-acid character AA."
  (sort (mapcar #'decode-aa-7bit
                (enum-encoded-aa (encode-aa-7bit aa)))
        #'char<=))

(defun enum-encoded-base (encoded-base)
  "Returns a list of all the unambiguous encoded bases represented by
ENCODED-BASE."
  (loop
     for b from 0 below (integer-length encoded-base)
     when (logbitp b encoded-base)
     collect (ash 1 b)))

(defun enum-encoded-codon (encoded-codon)
  "Returns a list of all the unambiguous encoded codons represented by
ENCODED-CODON."
  (if (null (rest encoded-codon))
      (mapcar #'list (enum-encoded-base (first encoded-codon)))
    (loop
       for x in (enum-encoded-base (first encoded-codon))
       nconc (mapcar (lambda (y)
                       (cons x y))
                     (enum-encoded-codon (rest encoded-codon))))))

(defun enum-encoded-aa (encoded-aa)
   "Returns a list of all the unambiguous encoded amino-acids
represented by ENCODED-AA."
  (cond ((= 2 encoded-aa)  ; aspartatic-acid/asparagine
         (list 4 14))
        ((= 10 encoded-aa) ; isoleucine/leucine
         (list 9 12))
        ((= 26 encoded-aa) ; glutamic-acid/glutamine
         (list 1 17))
        (t
         (list encoded-aa))))

(defun phred-quality (p)
  "Returns the Phred score of a base where P is the error
probability."
  (round (* -10 (/ (log p)
                   (log 10)))))

(defun phred-probability (q)
  "Returns the error probability of a base where Q is the Phred
score."
  (expt 10 (/ (- q) 10)))

(defun encode-phred-quality (q)
  "Returns the character encoding of the Phred quality score Q."
  (code-char (+ 33 (min 93 q))))

(defun decode-phred-quality (c)
  "Returns the Phred quality score encoded by the character C."
  (- (char-code c) 33))

(defun illumina-quality (p)
  "Returns the Illumina score of a base where P is the error
probability."
  (let ((denom (- 1 p)))
    (if (zerop denom)
        -20
      (round (* -10 (/ (log (/ p denom))
                       (log 10)))))))

(defun encode-illumina-quality (q)
  "Returns a character encoding of the Illumina quality score Q."
  (code-char (+ 64 q)))

(defun decode-illumina-quality (c)
  "Returns the Illumina quality score encoded by the character C."
  (- (char-code c) 64))

(defun illumina-to-phred-quality (q)
  "Returns the Phred quality score corresponding to the Illumina
quality score Q."
  (* 10 (/ (log (1+ (expt 10 (/ q 10))))
           (log 10))))
