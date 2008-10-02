;;;
;;; Copyright (C) 2007-2008, Keith James. All rights reserved.
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

(defun complement-dna (base)
  "Returns the complentary DNA base to that represented by the
character BASE."
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
character BASE."
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

(declaim (inline complement-dna-4bit))
(defun complement-dna-4bit (encoded-base)
  "Returns the complentary encoded DNA base to that represented by
ENCODED-BASE."
  (ecase encoded-base
    (#b0001 #b0100)
    (#b0010 #b1000)
    (#b0100 #b0001)
    (#b1000 #b0010)
    (#b1100 #b0011)
    (#b0011 #b1100)
    (#b1001 #b0110)
    (#b0110 #b1001)
    (#b1010 #b1010)
    (#b0101 #b0101)
    (#b1011 #b1110)
    (#b1101 #b0111)
    (#b0111 #b1101)
    (#b1110 #b1011)
    (#b1111 #b1111)
    (#b0000 #b0000)))

(declaim (inline complement-rna-4bit))
(defun complement-rna-4bit (encoded-base)
  "Returns the complentary encoded RNA base to that represented by
ENCODED-BASE."
  (ecase encoded-base
    (#b0001 #b0100)
    (#b0010 #b1000)
    (#b0100 #b0001)
    (#b1000 #b0010)
    (#b1100 #b0011)
    (#b0011 #b1100)
    (#b1001 #b0110)
    (#b0110 #b1001)
    (#b1010 #b1010)
    (#b0101 #b0101)
    (#b1011 #b1110)
    (#b1101 #b0111)
    (#b0111 #b1101)
    (#b1110 #b1011)
    (#b1111 #b1111)
    (#b0000 #b0000)))

(declaim (inline encode-dna-4bit))
(defun encode-dna-4bit (base)
  "Encodes DNA standard-char BASE as a 4-bit byte, representing T as
0001, C as 0010, A as 0100 and G as 1000. The first base is in the
most significant 4-bit byte and the last base is in the least
significant 4-bit byte. Ambiguous bases are represented by bitwise AND
combinations of these."
  (ecase base
    ((#\t #\T) #b0001)
    ((#\c #\C) #b0010)
    ((#\a #\A) #b0100)
    ((#\g #\G) #b1000)
    ((#\r #\R) #b1100)
    ((#\y #\Y) #b0011)
    ((#\k #\K) #b1001)
    ((#\m #\M) #b0110)
    ((#\s #\S) #b1010)
    ((#\w #\W) #b0101)
    ((#\b #\B) #b1011)
    ((#\d #\D) #b1101)
    ((#\h #\H) #b0111)
    ((#\v #\V) #b1110)
    ((#\n #\N) #b1111)
    (#\-       #b0000)))

(declaim (inline decode-dna-4bit))
(defun decode-dna-4bit (encoded-base)
  "Decodes 4-bit byte ENCODED-BASE, returning a lower case character."
  (ecase encoded-base
    (#b0001 #\t)
    (#b0010 #\c)
    (#b0100 #\a)
    (#b1000 #\g)
    (#b1100 #\r)
    (#b0011 #\y)
    (#b1001 #\k)
    (#b0110 #\m)
    (#b1010 #\s)
    (#b0101 #\w)
    (#b1011 #\b)
    (#b1101 #\d)
    (#b0111 #\h)
    (#b1110 #\v)
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
    ((#\u #\U) #b0001)
    ((#\c #\C) #b0010)
    ((#\a #\A) #b0100)
    ((#\g #\G) #b1000)
    ((#\r #\R) #b1100)
    ((#\y #\Y) #b0011)
    ((#\k #\K) #b1001)
    ((#\m #\M) #b0110)
    ((#\s #\S) #b1010)
    ((#\w #\W) #b0101)
    ((#\b #\B) #b1011)
    ((#\d #\D) #b1101)
    ((#\h #\H) #b0111)
    ((#\v #\V) #b1110)
    ((#\n #\N) #b1111)
    (#\-       #b0000)))

(declaim (inline decode-rna-4bit))
(defun decode-rna-4bit (encoded-base)
  "Decodes 4-bit byte ENCODED-BASE, returning a lower case character."
  (ecase encoded-base
    (#b0001 #\u)
    (#b0010 #\c)
    (#b0100 #\a)
    (#b1000 #\g)
    (#b1100 #\r)
    (#b0011 #\y)
    (#b1001 #\k)
    (#b0110 #\m)
    (#b1010 #\s)
    (#b0101 #\w)
    (#b1011 #\b)
    (#b1101 #\d)
    (#b0111 #\h)
    (#b1110 #\v)
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
  (code-char (+ 33 (if (< 93 q)
                       93
                     q))))

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
