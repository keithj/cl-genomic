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
(defconstant +insertx+ 1
  "Dynamic programming backtrace delete pointer.")
(defconstant +inserty+ 2
  "Dynamic programming backtrace insert pointer.")
(defconstant +match+ 3
  "Dynamic programming backtrace match pointer.")
(defconstant +none+ 0
  "Dynamic programming backtrace null pointer.")

(defmacro with-affine-gap-matrices ((score insertx inserty backtrace) (m n)
                                    &body body)
  `(let ((,score (make-array (list ,m ,n) :element-type 'single-float
                             :initial-element 0.0))
         (,insertx (make-array (list ,m ,n) :element-type 'single-float
                               :initial-element 0.0))
         (,inserty (make-array (list ,m ,n) :element-type 'single-float
                               :initial-element 0.0))
         (,backtrace (make-array (list ,m ,n) :element-type 'path-pointer
                                 :initial-element 0)))
    (declare (type (simple-array single-float (* *)) ,score ,insertx ,inserty)
             (type (simple-array path-pointer (* *)) ,backtrace))
    ,@body))

(defmacro define-affine-gap-dp (((cell prev-cell max-cell)
                                 (cell-score max-score)
                                 (score insertx inserty backtrace
                                        (gap-open gap-extend))
                                 subst-form &optional cell-exclusion-form)
                                &body body)
  (destructuring-bind ((row col) (prev-row prev-col) (max-row max-col))
      (list cell prev-cell max-cell)
    (with-gensyms (rows cols ix-score iy-score ss
                        ix1 ix2 ix3 iy1 iy2 iy3 s1 s2 s3 ms)
      `(let ((,rows (array-dimension ,score 0))
             (,cols (array-dimension ,score 1))
             (,max-score 0.0)
             (,max-row 0)
             (,max-col 0))
        (loop
           for ,row of-type fixnum from 1 below ,rows
           for ,prev-row of-type fixnum = (1- ,row)
           do
             (loop
                for ,col of-type fixnum from 1 below ,cols
                for ,prev-col of-type fixnum = (1- ,col)
                  ,@(when cell-exclusion-form
                          `(when ,cell-exclusion-form))
                do
                  (let ((,ss ,subst-form))
                    (let ((,ix1 (+ (aref ,score ,row ,prev-col) ,gap-open))
                          (,ix2 (+ (aref ,insertx ,row ,prev-col) ,gap-extend))
                          (,ix3 (+ (aref ,inserty ,row ,prev-col) ,gap-open))
                          (,iy1 (+ (aref ,score ,prev-row ,col) ,gap-open))
                          (,iy2 (+ (aref ,inserty ,prev-row ,col) ,gap-extend))
                          (,iy3 (+ (aref ,insertx ,prev-row ,col) ,gap-open))
                          (,s1 (+ (aref ,score ,prev-row ,prev-col) ,ss))
                          (,s2 (+ (aref ,insertx ,prev-row ,prev-col) ,ss))
                          (,s3 (+ (aref ,inserty ,prev-row ,prev-col) ,ss)))
                      (let ((,ix-score (if (= 1 ,col)
                                           ,ix1
                                         (max ,ix1 ,ix2 ,ix3)))
                            (,iy-score (if (= 1 ,row)
                                           ,iy1
                                         (max ,iy1 ,iy2 ,iy3)))
                            (,cell-score (max 0.0 ,s1 ,s2 ,s3)))
                        (setf (aref ,insertx ,row ,col) ,ix-score
                              (aref ,inserty ,row ,col) ,iy-score
                              (aref ,score ,row ,col) ,cell-score)
                        (let ((,ms (max ,cell-score ,ix-score ,iy-score)))
                          (setf (aref ,backtrace ,row ,col)
                                (cond ((= ,ms ,ix-score)
                                       +insertx+)
                                      ((= ,ms ,iy-score)
                                       +inserty+)
                                      (t
                                       +match+))))
                        (when (> ,cell-score ,max-score)
                          (setf ,max-score ,cell-score
                                ,max-row ,row
                                ,max-col ,col)))))))
        ,@body))))

(defmacro define-backtrace (((vecm vecn alm aln &key gap)
                             (score insertx inserty)
                             (btrace (start-row start-col) (row col)))
                            &body body)
  (with-gensyms (ptr)
    `(let ((,row ,start-row)
           (,col ,start-col))
      (declare (type (simple-array path-pointer (* *)) ,btrace))
      (loop
         while (and (plusp (aref ,btrace ,row ,col))
                (or (plusp (aref ,score ,row ,col))
                    (plusp (aref ,insertx ,row ,col))
                    (plusp (aref ,inserty ,row ,col))))
         do (let ((,ptr (aref ,btrace ,row ,col)))
            (cond ((= +match+ ,ptr)
                   (vector-push-extend (aref ,vecm (1- ,row)) ,alm)
                   (vector-push-extend (aref ,vecn (1- ,col)) ,aln)
                   (decf ,row)
                   (decf ,col))
                  ((= +insertx+ ,ptr)
                   (vector-push-extend (aref ,vecn (1- ,col)) ,aln)
                   (vector-push-extend ,gap ,alm)
                   (decf ,col))
                  ((= +inserty+ ,ptr)
                   (vector-push-extend (aref ,vecm (1- ,row)) ,alm)
                   (vector-push-extend ,gap ,aln)
                   (decf ,row))
                  (t
                   (error "Unknown backtrace pointer ~a" ,ptr)))))
      ,@body)))

;;; The alignment return values will be refined when we have a proper
;;; sequence alignment representation
(defmethod align-local ((seqm encoded-dna-sequence)
                        (seqn encoded-dna-sequence) subst-fn
                        &key (gap-open -5.0) (gap-extend -1.0)
                        (band-centre 0)
                        (band-width most-positive-fixnum) alignment)
  (multiple-value-bind (align-score align)
      (smith-waterman-gotoh-4bit (vector-of seqm) (vector-of seqn) subst-fn
                                 :gap-open gap-open :gap-extend gap-extend
                                 :band-centre band-centre
                                 :band-width band-width
                                 :alignment alignment)
    (values align-score align)))

;; (defmethod align-local ((seqm dna-quality-sequence)
;;                         (seqn dna-quality-sequence) subst-fn
;;                         &key (gap-open -10.0) (gap-extend -1.0)
;;                         (band-centre 0)
;;                         (band-width most-positive-fixnum) alignment)
;;   (declare (ignore band-centre band-width))  
;;   (smith-waterman-gotoh-qual (vector-of seqm) (vector-of seqn)
;;                              (quality-of seqm) (quality-of seqn) subst-fn
;;                              :gap-open gap-open :gap-extend gap-extend
;;                              :alignment alignment))

(defmethod align-local ((seqm encoded-aa-sequence)
                        (seqn encoded-aa-sequence) subst-fn
                        &key (gap-open -10.0) (gap-extend -1.0)
                        (band-centre 0)
                        (band-width most-positive-fixnum) alignment)
  (multiple-value-bind (align-score align)
      (smith-waterman-gotoh-7bit (vector-of seqm) (vector-of seqn) subst-fn
                                 :gap-open gap-open :gap-extend gap-extend
                                 :band-centre band-centre
                                 :band-width band-width
                                 :alignment alignment)
    (values align-score align)))

;; Modify finding shared kmers to use a substitution matrix
;; kmer seeded heuristic
(defmethod align-local-ksh ((seqm encoded-dna-sequence)
                            (seqn encoded-dna-sequence) subst-fn
                            &key (k 6) (gap-open -5.0) (gap-extend -1.0)
                            alignment)
  (let ((vecm (vector-of seqm))
        (vecn (vector-of seqn))
        (align-score 0.0)
        (align nil))
    (multiple-value-bind (bwidth bcentre)
        (multiple-value-call #'pairwise-band-width
          (find-common-kmers vecm vecn k))
      (when (plusp bwidth)
        (multiple-value-setq (align-score align)
          (smith-waterman-gotoh-4bit vecm vecn subst-fn
                                     :gap-open gap-open :gap-extend gap-extend
                                     :band-centre bcentre :band-width bwidth
                                     :alignment alignment))))
    (values align-score align)))

(defun smith-waterman-gotoh-4bit (vecm vecn subst-fn
                                  &key (gap-open -5.0)
                                  (gap-extend -1.0)
                                  (band-centre 0)
                                  (band-width most-positive-fixnum)
                                  alignment)
  (declare (optimize (speed 3) (safety 0)))
  (declare (type function subst-fn)
           (type single-float gap-open gap-extend)
           (type fixnum band-centre band-width)
           (type (simple-array (unsigned-byte 4) (*)) vecm vecn))
  (flet ((subn (x y) ; local fn to avoid boxing of returned floats
           (funcall subst-fn x y)))
    (let ((m (length vecm))
          (n (length vecn))
          (half-width (ceiling band-width 2)))
      (with-affine-gap-matrices
          (sc ix iy bt) ((1+ m) (1+ n))
        (define-affine-gap-dp
            (((row col) (prev-row prev-col) (max-row max-col))
             (cell-score max-score)
             (sc ix iy bt (gap-open gap-extend))
             (the single-float
               (subn (aref vecm prev-row) (aref vecn prev-col)))
             (let ((diag (- row col)))
               (< (- diag half-width) band-centre (+ diag half-width))))
          (values max-score
                  (when alignment
                    (multiple-value-bind (alm startm endm aln startn endn)
                        (dp-backtrace-4bit vecm vecn sc ix iy
                                           bt max-row max-col)
                      (make-instance 'alignment
                                     :intervals (list
                                                 (make-na-align-interval
                                                  alm startm endm)
                                                 (make-na-align-interval
                                                  aln startn endn)))))))))))

(defun smith-waterman-gotoh-7bit (vecm vecn subst-fn
                                  &key (gap-open -10.0)
                                  (gap-extend -1.0)
                                  (band-centre 0)
                                  (band-width most-positive-fixnum)
                                  alignment)
  (declare (optimize (speed 3) (safety 0)))
  (declare (type function subst-fn)
           (type single-float gap-open gap-extend)
           (type fixnum band-centre band-width)
           (type (simple-array (unsigned-byte 7) (*)) vecm vecn))
  (flet ((subn (x y) ; local fn to avoid boxing of returned floats
           (funcall subst-fn x y)))
    (let ((m (length vecm))
          (n (length vecn))
          (half-width (ceiling band-width 2)))
      (with-affine-gap-matrices
          (sc ix iy bt) ((1+ m) (1+ n))
        (define-affine-gap-dp
            (((row col) (prev-row prev-col) (max-row max-col))
             (cell-score max-score)
             (sc ix iy bt (gap-open gap-extend))
             (the single-float
               (subn (aref vecm prev-row) (aref vecn prev-col)))
             (let ((diag (- row col)))
               (< (- diag half-width) band-centre (+ diag half-width))))
          (values max-score
                  (when alignment
                    (multiple-value-bind (alm startm endm aln startn endn)
                        (dp-backtrace-7bit vecm vecn sc ix iy
                                           bt max-row max-col)
                      (make-instance 'alignment
                                     :intervals (list
                                                 (make-aa-align-interval
                                                  alm startm endm)
                                                 (make-aa-align-interval
                                                  aln startn endn)))))))))))

;; (defun smith-waterman-gotoh-qual (vecm vecn qualm qualn subst-fn
;;                                   &key (gap-open -10.0) (gap-extend -1.0)
;;                                   alignment)
;;   (flet ((subn (x y qx qy)
;;            (funcall subst-fn x y qx qy)))
;;     (let ((m (length vecm))
;;           (n (length vecn)))
;;       (with-affine-gap-matrices
;;           (mat del ins btr) ((1+ m) (1+ n))
;;         (define-affine-gap-dp
;;             ( ((row col) (prev-row prev-col) (max-row max-col))
;;               (cell-score del-score ins-score max-score)
;;               (mat del ins btr (gap-open gap-extend))
;;              (max 0.0
;;                   (+ (aref mat prev-row prev-col)
;;                      (subn (aref vecm prev-row) (aref vecn prev-col)
;;                            (aref qualm prev-row) (aref qualn prev-col)))))
;;           (values max-score
;;                   (when alignment
;;                     (dp-backtrace vecm vecn mat btr max-row max-col))))))))

(defun dp-backtrace-4bit (vecm vecn score insertx inserty
                          btrace start-row start-col)
  (declare (optimize (speed 3) (safety 0)))
  (declare (type (simple-array (unsigned-byte 4) (*)) vecm vecn)
           (type (simple-array single-float (* *)) score insertx inserty))
  (let ((alm (make-array 100 :element-type '(unsigned-byte 4)
                         :adjustable t :fill-pointer 0))
        (aln (make-array 100 :element-type '(unsigned-byte 4)
                             :adjustable t :fill-pointer 0)))
    (define-backtrace ((vecm vecn alm aln :gap 0)
                       (score insertx inserty)
                       (btrace (start-row start-col) (row col)))
      (values (nreverse alm) row col (nreverse aln) start-row start-col))))

(defun dp-backtrace-7bit (vecm vecn score insertx inserty
                          btrace start-row start-col)
  (declare (optimize (speed 3) (safety 0)))
  (declare (type (simple-array (unsigned-byte 7) (*)) vecm vecn)
           (type (simple-array single-float (* *)) score insertx inserty))
  (let ((alm (make-array 100 :element-type '(unsigned-byte 7)
                             :adjustable t :fill-pointer 0))
        (aln (make-array 100 :element-type '(unsigned-byte 7)
                             :adjustable t :fill-pointer 0)))
    (define-backtrace ((vecm vecn alm aln :gap 0)
                       (score insertx inserty)
                       (btrace (start-row start-col) (row col)))
      (values (nreverse alm) row col (nreverse aln) start-row start-col))))

(defun make-na-align-interval (vec lower upper
                               &optional (strand *forward-strand*))
  (make-instance 'na-alignment-interval
                 :lower lower :upper upper :strand strand
                 :aligned (make-dna vec)))

(defun make-aa-align-interval (vec lower upper)
  (make-instance 'aa-alignment-interval
                 :lower lower :upper upper :aligned (make-aa vec)))

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
                (align-local adapter
                             (make-dna (assocdr :residues fq))
                             #'simple-dna-subst)
              ;; (align-local-ksh adapter (make-dna (assocdr :residues fq))
              ;;                 #'simple-dna-subst :k 6)
              (when (> score 35.0)
                (incf total))
              (when (= 50000 count)
                 (return))
              (when (zerop (rem count 1000))
                (format t "~a ... ~a~%" count seqs))
              (incf count))
         finally (return total)))))

(defun read-fasta (filespec)
  (with-open-file (fs filespec
                   :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (next (make-seq-input (make-line-input-stream fs) :fasta
                          :alphabet :dna))))
