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

(defmacro with-affine-gap-matrices ((score delete insert backtrace) (m n)
                                    &body body)
  `(let ((,score (make-array (list ,m ,n) :element-type 'single-float
                             :initial-element 0.0))
         (,delete (make-array (list ,m ,n) :element-type 'single-float
                              :initial-element 0.0))
         (,insert (make-array (list ,m ,n) :element-type 'single-float
                              :initial-element 0.0))
         (,backtrace (make-array (list ,m ,n) :element-type 'path-pointer
                                 :initial-element +none+)))
     (declare (type (simple-array single-float (* *))
                    ,score ,delete ,insert)
              (type (simple-array path-pointer (* *)) ,backtrace))
     ,@body))

(defmacro define-affine-gap-dp (((cell prev-cell max-cell)
                                 (cell-score del-score ins-score max-score)
                                 (score delete insert backtrace
                                        (gap-open gap-extend))
                                 cell-score-form &optional cell-exclusion-form)
                                &body body)
  (destructuring-bind ((row col) (prev-row prev-col) (max-row max-col))
      (list cell prev-cell max-cell)
    (with-gensyms (rows cols)
      `(let ((,rows (array-dimension ,score 0))
             (,cols (array-dimension ,score 1))
             (,max-score 0.0)
             (,max-row 0)
             (,max-col 0))
         (loop
            for ,row of-type fixnum from 1 below ,rows
            for ,prev-row of-type fixnum = (1- ,row)
            do (loop
                  for ,col of-type fixnum from 1 below ,cols
                  for ,prev-col of-type fixnum = (1- ,col)
                    ,@(when cell-exclusion-form
                        `(when ,cell-exclusion-form))
                  do (let ((,del-score
                            (if (= 1 ,row)
                                (+ (aref ,score ,prev-row ,col)
                                   ,gap-open)
                              (max (+ (aref ,score ,prev-row ,col)
                                      ,gap-open)
                                   (+ (aref ,delete ,prev-row ,col)
                                      ,gap-extend))))
                           (,ins-score
                            (if (= 1 ,col)
                                (+ (aref ,score ,row ,prev-col)
                                   ,gap-open)
                              (max (+ (aref ,score ,row ,prev-col)
                                      ,gap-open)
                                   (+ (aref ,insert ,row ,prev-col)
                                      ,gap-extend)))))
                       (setf (aref ,delete ,row ,col) ,del-score
                             (aref ,insert ,row ,col) ,ins-score)
                       (let ((,cell-score ,cell-score-form))
                         (setf (aref ,score ,row ,col) ,cell-score
                               (aref ,backtrace ,row ,col)
                               (cond ((= ,cell-score ,del-score)
                                      +delete+)
                                     ((= ,cell-score ,ins-score)
                                      +insert+)
                                     (t
                                      +match+)))
                         (when (> ,cell-score ,max-score)
                           (setf ,max-score ,cell-score
                                 ,max-row ,row
                                 ,max-col ,col))))))
         ,@body))))

;; Modify finding shared kmers to use a substitution matrix

;; kmer seeded heuristic
(defmethod align-local-ksh ((seqm encoded-dna-sequence)
                            (seqn encoded-dna-sequence) subst-fn
                            &key (k 6) (gap-open -10.0) (gap-extend -1.0)
                            alignment)
  (let ((vecm (vector-of seqm))
        (vecn (vector-of seqn))
        (align-score 0.0)
        (align-obj nil))
    (multiple-value-bind (bwidth bcentre)
        (multiple-value-call #'pairwise-band-width
          (find-common-kmers vecm vecn k))
      (when (plusp bwidth)
        (multiple-value-setq (align-score align-obj)
          (smith-waterman-gotoh vecm vecn subst-fn
                                :gap-open gap-open :gap-extend gap-extend
                                :band-centre bcentre :band-width bwidth
                                :alignment alignment))))
    (values align-score align-obj)))

(defmethod align-local-ksh ((seqm vector) (seqn vector) subst-fn
                            &key (k 6) (gap-open -10.0) (gap-extend -1.0)
                            alignment)
  (let ((vecm (ensure-encoded-4bit seqm #'encode-dna-4bit))
        (vecn (ensure-encoded-4bit seqn #'encode-dna-4bit))
        (align-score 0.0)
        (align-obj nil))
    (multiple-value-bind (bwidth bcentre)
        (multiple-value-call #'pairwise-band-width
          (find-common-kmers vecm vecn k))
      (when (plusp bwidth)
        (multiple-value-setq (align-score align-obj)
          (smith-waterman-gotoh vecm vecn subst-fn
                                :gap-open gap-open :gap-extend gap-extend
                                :band-centre bcentre :band-width bwidth
                                :alignment alignment))))
    (values align-score align-obj)))

(defmethod align-local ((seqm encoded-dna-sequence)
                        (seqn encoded-dna-sequence) subst-fn
                        &key (gap-open -10.0) (gap-extend -1.0)
                        (band-centre 0)
                        (band-width most-positive-fixnum) alignment)
  (smith-waterman-gotoh (vector-of seqm) (vector-of seqn) subst-fn
                        :gap-open gap-open :gap-extend gap-extend
                        :band-centre band-centre :band-width band-width
                        :alignment alignment))

(defmethod align-local ((seqm vector) (seqn vector) subst-fn
                        &key (gap-open -10.0) (gap-extend -1.0)
                        (band-centre 0)
                        (band-width most-positive-fixnum) alignment)
  (let ((vecm (ensure-encoded-4bit seqm #'encode-dna-4bit))
        (vecn (ensure-encoded-4bit seqn #'encode-dna-4bit)))
    (smith-waterman-gotoh vecm vecn subst-fn
                          :gap-open gap-open :gap-extend gap-extend
                          :band-centre band-centre :band-width band-width
                          :alignment alignment)))

(defmethod align-local ((seqm dna-quality-sequence)
                        (seqn dna-quality-sequence) subst-fn
                        &key (gap-open -10.0) (gap-extend -1.0)
                        (band-centre 0)
                        (band-width most-positive-fixnum) alignment)
  (declare (ignore band-centre band-width))  
  (smith-waterman-gotoh-qual (vector-of seqm) (vector-of seqn)
                             (quality-of seqm) (quality-of seqn) subst-fn
                             :gap-open gap-open :gap-extend gap-extend
                             :alignment alignment))

(defun smith-waterman-gotoh (vecm vecn subst-fn
                             &key (gap-open -10.0) (gap-extend -1.0)
                             (band-centre 0) (band-width most-positive-fixnum)
                             alignment)
  "Performs a Smith Waterman alignment of VECM against VECN using
Gotoh's improvement. The affine gap scoring is expressed as penalties,
that is to say GAP-OPEN and GAP-EXTEND should be negative values. The
alignment may be banded to prune the search space.

Arguments:

- vecm \(vector\): A vector to be aligned.
- vecn \(vector\): A vector to be aligned.
- subst-fn \(function\): A substitution function that accepts two
vector elements as arguments and returns a single-float substitution
score.

Key:

- gap-open \(single-float\): The gap opening penalty.
- gap-extend \(single-float\): The gap extension penalty.
- band-centre \(fixnum\): The band centre for banded searches. This
defaults to 0, the main diagonal. The desired band may be calculated
by subtracting a j coordinate from its corresponding i coordinate.
- band-width \(fixnum\): The band width for banded searches. This
defaults to most-positive-fixnum so that the search space is not
pruned.
- alignment \(boolean\): Flag to indicate whether an alignment should
be returned.

Returns:

- The single-float score of the alignment.
- An alignment object, or NIL."
  (declare (optimize (speed 3) (safety 0)))
  (declare (type function subst-fn)
           (type single-float gap-open gap-extend)
           (type fixnum band-centre band-width)
           (type (simple-array (unsigned-byte 4)) vecm vecn))
  (flet ((subn (x y) ; local fn to avoid boxing of returned floats
           (funcall subst-fn x y)))
    (let ((m (length vecm))
          (n (length vecn))
          (half-width (ceiling band-width 2)))
      (with-affine-gap-matrices
          (mat del ins btr) ((1+ m) (1+ n))
        (define-affine-gap-dp
            (((row col) (prev-row prev-col) (max-row max-col))
             (cell-score del-score ins-score max-score)
             (mat del ins btr (gap-open gap-extend))
             (max 0.0 ; to keep the search local
                  (the single-float
                    (+ (aref mat prev-row prev-col)
                       (the single-float
                         (subn (aref vecm prev-row) (aref vecn prev-col))))))
             (let ((diag (- row col)))
               (< (- diag half-width) band-centre (+ diag half-width))))
          (values max-score
                  (when alignment
                    (dp-backtrace vecm vecn mat btr max-row max-col))))))))

(defun smith-waterman-gotoh-qual (vecm vecn qualm qualn subst-fn
                                  &key (gap-open -10.0) (gap-extend -1.0)
                                  alignment)
  (flet ((subn (x y qx qy)
           (funcall subst-fn x y qx qy)))
    (let ((m (length vecm))
          (n (length vecn)))
      (with-affine-gap-matrices
          (mat del ins btr) ((1+ m) (1+ n))
        (define-affine-gap-dp
            ( ((row col) (prev-row prev-col) (max-row max-col))
              (cell-score del-score ins-score max-score)
              (mat del ins btr (gap-open gap-extend))
             (max 0.0
                  (+ (aref mat prev-row prev-col)
                     (subn (aref vecm prev-row) (aref vecn prev-col)
                           (aref qualm prev-row) (aref qualn prev-col)))))
          (values max-score
                  (when alignment
                    (dp-backtrace vecm vecn mat btr max-row max-col))))))))

(defun dp-backtrace (vecm vecn scomat btmat bti btj)
  (declare (optimize (speed 3) (safety 0)))
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
    (list i (make-instance 'encoded-dna-sequence :vector (nreverse am)) bti
          j (make-instance 'encoded-dna-sequence :vector (nreverse an)) btj)))

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
  (declare (optimize (speed 3) (safety 0)))
  (declare (type (simple-array (unsigned-byte 4) (*)) vec)
           (type fixnum k))
  (let ((end (- (length vec) k))
        (table (make-hash-table :size size :test #'equalp)))
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
  (declare (optimize (speed 3) (safety 0)))
  (declare (type (simple-array (unsigned-byte 4) (*)) vecm vecn)
           (type fixnum k))
  (let ((end (- (length vecm) k))
        (kmern (make-kmer-table vecn k size))
        (matched (make-hash-table :size size :test #'equalp)))
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
  (declare (optimize (speed 3) (safety 0)))
  (declare (type list mcoords ncoords))
  (let ((k 0)
        (c 0))
    (declare (type fixnum k c))
    (when (and mcoords ncoords)
      (let ((mind 0)
            (maxd 0))
        (declare (type fixnum mind maxd))
        (loop
           for m in mcoords
           for n in ncoords
           do (loop
                 for i of-type fixnum in m
                 do (loop
                       for j of-type fixnum in n
                       for d of-type fixnum = (- i j)
                       do (setf mind (min mind d)
                                maxd (max maxd d)))))
        ;; FIXME -- is this value for k correct? When there are no
        ;; kmers k is 2, which is why we test of that condition
        ;; above. It would be nice to omit that special case.
        (setf k (+ 2 (- maxd mind))
              c (round (+ (/ (the fixnum (- maxd mind)) 2.0)
                          mind)))))
    (values k c)))

(defun print-banding (m n band-centre band-width a b)
  (let ((x (make-array (list (1+ m) (1+ n)) :initial-element 0))
        (half-width (ceiling band-width 2)))
    (loop
       for i from 1 to m
       do (loop
             for j from 1 to n
             for diag = (- i j)
             do (progn
                  (when (< (- diag half-width)
                           band-centre
                           (+ diag half-width))
                    (setf (aref x i j) 1))
                  (when (= band-centre diag)
                    (setf (aref x i j) 8)))))
    (mapc #'(lambda (z w)
              (setf (aref x
                          (first z)
                          (first w)) 2)) a b)
    (princ x)
    (terpri)))

(defun test-align (fastq-filespec)
  (with-open-file (in fastq-filespec
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let ((adapter (make-dna "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"))
          (gen  (make-seq-input (make-line-input-stream in) :fastq
                                :alphabet :dna :metric :phred
                                :parser (make-instance 'raw-sequence-parser))))
      (loop
         with total = 0
         with count = 0
         as fq = (next gen)
         while fq
         do (multiple-value-bind (score seqs)
                ;; (align-local (to-string adapter)
                ;;              (assocdr :residues fq) #'simple-dna-subst)
              (align-local-ksh (to-string adapter) (assocdr :residues fq)
                               #'simple-dna-subst :k 6)
              (when (> score 35.0)
                (incf total))
              (when (= 50000 count)
                (return))
              (when (zerop (rem count 1000))
                (format t "~a ...~%" count))
              (incf count))
         finally (return total)))))

(defun read-fasta (filespec)
  (with-open-file (fs filespec
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (next (make-seq-input (make-line-input-stream fs) :fasta
                          :alphabet :dna))))


