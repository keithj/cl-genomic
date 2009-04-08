
(in-package :cl-user)

(defparameter *tpm-thresholds*
  (list 1 5 10 50 100 500 1000 5000 10000))

(defun tag-abundance-fn (x)
  "Empirical tag abundance function calculated from human cell line
MPSS by Jongeneel et al. 2003."
  (exp (/ (* -2 (expt x 2)) (1+ x))))

(defun prop-genes-per-threshold (thresholds)
  "Returns the proportion of the total number of genes expected per
tpm threshold, given a list of THRESHOLDS."
  (let ((lts (mapcar #'(lambda (x)
                         (log x 10)) thresholds)))
    (loop for cons on lts
       if (second cons)
       collect (- (tag-abundance-fn (first cons))
                  (tag-abundance-fn (second cons)))
       else
       collect (tag-abundance-fn (first cons)))))

(defun cumul-pc-genes-above (thresholds)
  "Returns a list of percentage values that reflect how many genes
have at least that abundance."
  (maplist #'(lambda (x)
               (reduce #'+ x))
           (mapcar #'(lambda (x)
                       (* 100 x))
                   (prop-genes-per-threshold thresholds))))

(defun genes-per-threshold (ngenes thresholds)
  "Returns the number of genes per tpm threshold, given a list of
THRESHOLDS and the total number of genes NGENES."
  (let ((counts (mapcar #'(lambda (x)
                            (round (* ngenes x)))
                        (prop-genes-per-threshold thresholds))))
    ;; The final count is made by subtraction from NGENES to ensure
    ;; the total is equal to NGENES
    (nconc (butlast counts)
           (list (- ngenes (reduce #'+ (butlast counts)))))))

(defun transcripts-per-threshold (ngenes thresholds)
  "Returns the number of transcripts per tpm threshold, given a list
of THRESHOLDS and the total number of genes NGENES."
  (mapcar #'* (genes-per-threshold ngenes thresholds) thresholds))

;; Result is 1395316.
;; (reduce #'+ (transcripts-per-threshold 30000 *tpm-thresholds*))

(defun rand-informative-p (prob)
  (<= (random 1.0) prob))

(defun make-sampler (n m)
  "Returns a closure that will randomly sample N integers in the range
1 to M, without replacement."
  (when (> n m)
    (error "Invalid n ~a (greater than m ~a)." n m))
  (let ((sampled (make-hash-table)))
    (lambda ()
      (if (= n (hash-table-count sampled))
          nil
        (let ((sample (loop
                         for s = (random m) then (random m)
                         while (gethash s sampled)
                         finally (return s))))
          (setf (gethash sample sampled) t)
          sample)))))

(defun sample-ests (nests ngenes thresholds &optional (iprob 1.0))
  (let ((genes-per-thresh (genes-per-threshold ngenes thresholds))
        (gene-id 0)
        (clone-id 0)
        (clones (make-hash-table)))
    (mapc #'(lambda (gene-count transcript-count)
              (loop for i of-type fixnum from 1 to gene-count
                 do (progn
                      (incf gene-id)
                      (loop for j of-type fixnum from 1 to transcript-count
                         do (progn
                              (incf clone-id)
                              (setf (gethash clone-id clones)
                                    (pairlis
                                     '(:category :gene-id :informative)
                                     (list transcript-count
                                           gene-id
                                           (rand-informative-p iprob)))))))))
          genes-per-thresh thresholds)
    (let ((sampler (make-sampler nests clone-id)))
      (loop
         for cid = (funcall sampler) then (funcall sampler)
         with clone
         while cid
         do (setf clone (gethash cid clones))
         when (cdr (assoc :informative clone))
         collect (cdr (assoc :gene-id clone)) into informative
         finally (return (remove-duplicates informative))))))


;; For an experiment, calulate the intensity values that define the
;; transcript abundances in categories above, assuming the
;; distribution given in TAG-ABUNDANCE-FN. Then take those thresholds
;; and apply them to to the intensities of differentially expressed
;; genes only, recording what fraction of them are in each
;; category. Finally, caculate the ratio of the two proportions for
;; each category.
