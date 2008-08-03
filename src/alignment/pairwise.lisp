;;;
;;; Copyright (C) 2008 Keith James. All rights reserved.
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

(deftype path-pointer ()
  "Dynamic programming backtrace pointer."
  '(unsigned-byte 2))
(defconstant +delete+ 1
  "Dynamic programming backtrace delete pointer.")
(defconstant +insert+ 2
  "Dynamic programming backtrace insert pointer.")
(defconstant +match+ 3
  "Dynamic programming backtrace match pointer.")
(defconstant +none+ 0
  "Dynamic programming backtrace null pointer.")


;; Modify finding shared kmers to use a substitution matrix

(defmethod align-local-affine ((seqm encoded-dna-sequence)
                               (seqn encoded-dna-sequence)
                               (submat subst-matrix)
                               &key (gap-open -10.0) (gap-extend -1.0)
                               (band-centre 0)
                               (band-width most-positive-fixnum)
                               alignment)
  (let ((vecm (vector-of seqm))
        (vecn (vector-of seqn)))
    (multiple-value-bind (max-score scomat btmat bti btj)
        (smith-waterman-gotoh vecm vecn
                              (matrix-of submat) (index-of submat)
                              :gap-open gap-open
                              :gap-extend gap-extend
                              :band-centre band-centre
                              :band-width band-width)
      (values max-score
              (when alignment
                (dp-backtrace vecm vecn scomat btmat bti btj))))))

(defun smith-waterman-gotoh (vecm vecn submat subidx &key
                             (gap-open -10.0) (gap-extend -1.0)
                             (band-centre 0) (band-width most-positive-fixnum))
  "Performs a Smith Waterman alignment of VECM against VECN using
Gotoh's improvement. Affine gap scoring is used, expressed as
penalties, that is to say GAP-OPEN and GAP-EXTEND should be negative
values. The alignment may be banded to prune the search space.

Arguments:

- vecm \(vector\): A vector to be aligned.
- vecn \(vector\): A vector to be aligned.
- submat \(matrix\): A substitution matrix.
- subidx \(vector\): A vector index used to locate rows and columns in
the substitution matrix. This index contains all the tokens in the
alphabet of the vectors being aligned. The index of the token in the
index corresponds to its row and column indices in the substitution
matrix.

Key:

- gap-open \(single-float\): The gap opening penalty.
- gap-extend \(single-float\): The gap extension penalty.

- band-centre \(fixnum\): The band centre for banded searches. This
defaults to 0, the main diagonal. The desired band may be calculated
by subtracting a j coordinate from its corresponding i coordinate.
- band-width \(fixnum\): The band width for banded searches. This
defaults to most-positive-fixnum so that the search space is not
pruned.

Returns:

- The single-float maximum score from the matrix.
- The single-float score matrix.
- The backtrace matrix.
- The i coordinate of the backtrace starting point.
- The j coordinate of the backtrace starting point."
  (declare (optimize (speed 3) (debug 0) (safety 1)))
  (declare (type single-float gap-open gap-extend)
           (type fixnum band-centre band-width)
           (type (simple-array single-float (* *)) submat)
           (type (simple-array (unsigned-byte 4)) vecm vecn subidx))
  (flet ((subn (x y) ; local fn to avoid boxing of returned floats
           (let ((i (position x subidx :test #'=))
                 (j (position y subidx :test #'=)))
             (aref submat i j))))
    (let* ((m (length vecm))
           (n (length vecn))
           (half-width (ceiling band-width 2))
           (dim (list (1+ m) (1+ n))))
      (let ((mat (make-array dim :element-type 'single-float
                             :initial-element 0.0))
            (del (make-array dim :element-type 'single-float
                             :initial-element 0.0))
            (ins (make-array dim :element-type 'single-float
                             :initial-element 0.0))
            (btr (make-array dim :element-type 'path-pointer
                             :initial-element +none+))
            (max-score 0.0)
            (imax 0)
            (jmax 0))
        (declare (type (simple-array single-float (* *)) mat del ins)
                 (type (simple-array path-pointer (* *)) btr))
        (loop
           for i of-type fixnum from 1 to m
           for idec of-type fixnum = (1- i)
           do (loop
                 for j of-type fixnum from 1 to n
                 for jdec of-type fixnum = (1- j)
                 for diag of-type fixnum = (- i j) ; the diagonal of this cell
                 when (< (- diag half-width) band-centre (+ diag half-width))
                 do (let ((dscore (if (= 1 i)
                                      (+ (aref mat idec j) gap-open)
                                    (max (+ (aref mat idec j) gap-open)
                                         (+ (aref del idec j) gap-extend))))
                          (iscore (if (= 1 j)
                                      (+ (aref mat i jdec) gap-open)
                                    (max (+ (aref mat i jdec) gap-open)
                                         (+ (aref ins i jdec) gap-extend)))))
                      (setf (aref del i j) dscore
                            (aref ins i j) iscore)
                      (let ((score (max 0.0 ; to keep the search local
                                        (+ (aref mat idec jdec)
                                           (subn (aref vecm idec)
                                                 (aref vecn jdec)))
                                        dscore
                                        iscore)))
                        (setf (aref mat i j) score
                              (aref btr i j) (cond ((= score dscore)
                                                    +delete+)
                                                   ((= score iscore)
                                                    +insert+)
                                                   (t
                                                    +match+)))
                        (when (> score max-score)
                          (setf max-score score
                                imax i
                                jmax j))))))
        (values max-score mat btr imax jmax)))))

(defun dp-backtrace (vecm vecn scomat btmat bti btj)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type (simple-array (unsigned-byte 4) (*)) vecm vecn)
           (type (simple-array single-float (* *)) scomat)
           (type (simple-array path-pointer (* *)) btmat))
  (let ((am (make-array 100 :element-type '(unsigned-byte 4)
                        :adjustable t :fill-pointer 0))
        (an (make-array 100 :element-type '(unsigned-byte 4)
                        :adjustable t :fill-pointer 0))
        (i bti)
        (j btj))
    (loop
       while (plusp (aref scomat i j))
       do (let ((ptr (aref btmat i j))
                (gap (encode-dna-4bit #\-)))
            (cond ((= +match+ ptr)
                   (vector-push-extend (aref vecm (1- i)) am)
                   (vector-push-extend (aref vecn (1- j)) an)
                   (decf i)
                   (decf j))
                  ((= +delete+ ptr)
                   (vector-push-extend (aref vecm (1- i)) am)
                   (vector-push-extend gap an)
                   (decf i))
                  ((= +insert+ ptr)
                   (vector-push-extend (aref vecn (1- j)) an)
                   (vector-push-extend gap am)
                   (decf j))
                  (t
                   (error "Unknown backtrace pointer ~a" ptr)))))
    (list (make-instance 'encoded-dna-sequence :vector am)
          (make-instance 'encoded-dna-sequence :vector an))))

(defun make-kmer-table (vec k &optional (size 16))
  "Creates a hash-table of the kmers of length K in vector VEC.

Arguments:
- vec (vector): A vector.
- k (fixnum): The kmer length.

Optional:
- size (fixnum): The initial hash-table size.

Returns:
A EQUAL hash-table of the kmers of length K in vector VEC. The hash
keys are the kmers, while the hash values are lists of the coordinates
at which those kmers occur."
  (declare (optimize (speed 3) (safety 1)))
  (declare (type simple-string vec)
           (type fixnum k))
  (let ((end (- (length vec) k))
        (table (make-hash-table :size size :test #'equal)))
    (loop
       for i from 0 to end
       do (push i (gethash (subseq vec i (+ i k)) table)))
    table))

(defun find-common-kmers (vecm vecn k &optional (size 16))
  "Finds the kmers of length K that occur in both VECM and VECN using
hash-tables.

Arguments:
- vecm (vector): A vector.
- vecn (vector): A vector.
- k (fixnum): The kmer length.

Optional:
- size (fixnum): The initial hash-table size.

Returns:
- A list of lists of fixnum start coordinates of kmers in VECM.
- A list of lists of fixnum start coordinates of kmers in VECN.

The returned lists each contain the a number of elements equal to the
number of unique, shared kmers found. The coordinate lists at the nth
position in each list refer to the same kmer."
  (declare (optimize (speed 3) (safety 1)))
  (declare (type simple-string vecm vecn)
           (type fixnum k))
  (let ((end (- (length vecm) k))
        (kmern (make-kmer-table vecn k size))
        (matched (make-hash-table :size size :test #'equal)))
    (loop
       for i from 0 to end
       as kmer = (subseq vecm i (+ i k))
       as j = (gethash kmer kmern)
       do (when j
            (push i (gethash kmer matched)))
       finally (return
                 (loop
                    for kmer being the hash-keys in matched
                    using (hash-value coords)
                    collect coords into mcoords
                    collect (gethash kmer kmern) into ncoords
                    finally (return (values mcoords ncoords)))))))

(defun pairwise-band-width (mcoords ncoords)
  "Calculates the matrix diagonal and minimum band width required to
prune a dynamic programming search to include all the kmer diagonals
described by the coordinates MCOORDS and NCOORDS.

Arguments:

- mcoords \(list list\): A list of lists of fixnum start coordinates
of kmers in sequence m.
- ncoords \(list list\): A list of lists of fixnum start coordinates,
the same length as mcoords, of the same kmers in sequence n.

The nth element in each argument list indicates the start coordinates
of the same kmer.

Returns:

- A fixnum matrix band width that contains all diagonals described by
the arguments.
- A fixnum diagonal about which the band is centred."
  (let ((mind 0)
        (maxd 0))
    (loop
       for m in mcoords
       for n in ncoords
       do (loop
             for i in m
             do (loop
                   for j in n
                   for d = (- i j)
                   do (setf mind (min mind d)
                            maxd (max maxd d)))))
    (let ((k (+ 2 (- maxd mind)))
          (c (round (+ (/ (- maxd mind) 2) mind))))
      (values k c))))


;; (defun test-align (fastq-filespec)
;;   (with-open-file (in fastq-filespec
;;                    :direction :input
;;                    :element-type 'base-char
;;                    :external-format :ascii)
;;     (let ((adapter "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG")
;;           (gen  (make-seq-input (make-line-input-stream in) :fastq
;;                                 :alphabet :dna :metric :phred
;;                                 :parser (make-instance
;;                                          'raw-sequence-parser))))
;;       (loop
;;          with total = 0
;;          with count = 0
;;          as fq = (next gen)
;;          while fq
;;          do (multiple-value-bind (score seqs)
;;                 (sw-affine3 adapter (assocdr :residues fq)
;;                             *simple-dna*
;;                             *simple-dna-index*
;;                             :gap-open -10.0  :gap-extend -1.0
;;                             :bandw 2)
;;               (when (> score 29)
;;                 (incf total)
;; ;;                 (princ score)
;; ;;                 (terpri)
;; ;;                 (write-line (first seqs))
;; ;;                 (write-line (second seqs))
;; ;;                 (terpri)
;;                 )
;;               (when (= 10000000 count)
;;                 (return))
;;               (when (zerop (rem count 100000))
;;                 (format t "~a ...~%" count))
;;               (incf count))
;;          finally (return total)))))

;; Find which area of the matrix contains all the interesting seeds
;; and limit the SW band to this area.
