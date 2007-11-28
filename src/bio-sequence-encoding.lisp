
(in-package :bio-sequence)

(declaim (inline encode-dna-2bit))
(defun encode-dna-2bit (base)
  "Encodes DNA standard-char BASE as a two-bit byte, representing T as
00, C as 01, A as 10 G and 11. The first base is in the most
significant two-bit byte and the last base is in the least significant
two-bit byte."
  (ecase base
    (#\t #b00)
    (#\c #b01)
    (#\a #b10)
    (#\g #b11)))

(declaim (inline decode-dna-2bit))
(defun decode-dna-2bit (encoded-residue)
  (ecase encoded-residue
    (#b00 #\t)
    (#b01 #\c)
    (#b10 #\a)
    (#b11 #\g)))

(declaim (inline encode-dna-4bit))
(defun encode-dna-4bit (base)
  "Encodes DNA standard-char BASE as a four-bit byte, representing T
as 0001, C as 0010, A as 0100 and G as 1000.  The first base is in the
most significant four-bit byte and the last base is in the least
significant four-bit byte. Ambiguous bases are represented by bitwise
AND combinations of these."
  (ecase base
    (#\t #b0001)
    (#\c #b0010)
    (#\a #b0100)
    (#\g #b1000)
    (#\r #b1100)
    (#\y #b0011)
    (#\k #b1001)
    (#\m #b0110)
    (#\s #b1010)
    (#\w #b0101)
    (#\b #b1011)
    (#\d #b1101)
    (#\h #b0111)
    (#\v #b1110)
    (#\n #b1111)))

(declaim (inline decode-dna-4bit))
(defun decode-dna-4bit (encoded-residue)
  (ecase encoded-residue
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
    (#b1111 #\n)))

(declaim (inline encode-rna-2bit))
(defun encode-rna-2bit (base)
  "Encodes RNA standard-char BASE as a two-bit byte, representing U as
00, C as 01, A as 10 G and 11. The first base is in the most
significant two-bit byte and the last base is in the least significant
two-bit byte."
  (ecase base
    (#\u #b00)
    (#\c #b01)
    (#\a #b10)
    (#\g #b11)))

(declaim (inline decode-rna-2bit))
(defun decode-rna-2bit (encoded-residue)
  (ecase encoded-residue
    (#b00 #\u)
    (#b01 #\c)
    (#b10 #\a)
    (#b11 #\g)))

(declaim (inline encode-rna-4bit))
(defun encode-rna-4bit (base)
  "Encodes RNA standard-char BASE as a four-bit byte, representing U
as 0001, C as 0010, A as 0100 and G as 1000.  The first base is in the
most significant four-bit byte and the last base is in the least
significant four-bit byte. Ambiguous bases are represented by bitwise
AND combinations of these."
  (ecase base
    (#\u #b0001)
    (#\c #b0010)
    (#\a #b0100)
    (#\g #b1000)
    (#\r #b1100)
    (#\y #b0011)
    (#\k #b1001)
    (#\m #b0110)
    (#\s #b1010)
    (#\w #b0101)
    (#\b #b1011)
    (#\d #b1101)
    (#\h #b0111)
    (#\v #b1110)
    (#\n #b1111)))

(declaim (inline decode-rna-4bit))
(defun decode-rna-4bit (encoded-residue)
  (ecase encoded-residue
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
    (#b1111 #\n)))

;; (declaim (inline encode-aa-8bit))
;; (defun encode-aa-8bit (aa)
;;   (ccase aa
;;     (#\A 0)
;;     (#\B 1)
;;     (#\C 2)
;;     (#\D 3)
;;     (#\E 4)
;;     (#\F 5)
;;     (#\G 6)
;;     (#\H 7)
;;     (#\I 8)
;;     (#\J 9)
;;     (#\K 10)
;;     (#\L 11)
;;     (#\M 12)
;;     (#\N 13)
;;     (#\O 14)
;;     (#\P 15)
;;     (#\Q 16)
;;     (#\R 17)
;;     (#\S 18)
;;     (#\T 19)
;;     (#\U 20)
;;     (#\V 21)
;;     (#\W 22)
;;     (#\X 23)
;;     (#\Y 24)
;;     (#\Z 25)))

(defun phred-quality (p)
  "Returns the Phred score of a base where P is the error
probability."
  (round (* -10 (/ (log p)
                   (log 10)))))

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
