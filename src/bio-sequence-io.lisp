
(in-package :bio-sequence)

(defparameter *seq-line-buffer-size* 512)

(define-condition malformed-record-error (error)
  ((text  :initform nil
          :initarg :text
          :reader text-of
          :documentation "Error message text."))
  (:report (lambda (condition stream)
             (format stream "A malformed record was encountered~@[: ~a~]"
                     (text-of condition)))))

(defun make-dna-sexp (name token-seq &optional description)
  "Returns a standard DNA sexp, given a sequence NAME, a vector of
residue tokens TOKEN-SEQ and a DESCRIPTION string."
  (cons :bio-seq (pairlis '(:alphabet :name :description
                            :ambiguity :token-seq)
                          (list :dna name description
                                (if (or (position (char-code #\N)
                                                  token-seq)
                                        (position (char-code #\n)
                                                  token-seq))
                                    :iupac
                                  nil)
                                token-seq))))

(defun make-quality-sexp (dna-sexp quality)
  "Returns a standard quality sexp, given a standard DNA-SEXP and a
quality vector."
  ;; FIXME -- check that token-seq and quality strings are the same
  ;; length
  (cons :quality-seq (list dna-sexp
                           (cons :quality quality))))

(defun make-removing-callback (callback predicate)
  (lambda (sexp &rest callback-args)
    (if (funcall predicate sexp)
        nil
      (apply callback sexp callback-args))))

(defun ambiguous-seq-p (sexp)
  (let ((content (cdr sexp)))
    (eql :iupac (cdr (assoc :ambiguity content)))))


(defun make-chunk-pname (file-pname chunk-count)
  (make-pathname
   :directory (pathname-directory file-pname)
   :name (concatenate 'string (pathname-name file-pname) "."
                      (princ-to-string chunk-count))
   :type (pathname-type file-pname)))


(defun benchmark (filename)
  (with-open-file (stream filename :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((line-buffer (make-line-buffer stream)))
      (do ((fq (read-fastq line-buffer)
               (read-fastq line-buffer))
           (old-fq nil fq)
           (counter 0 (1+ counter)))
          ((null fq) (list counter
                           old-fq))))))

(defun benchmark-fasta1 (filename)
  (with-open-file (stream filename :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((line-buffer (make-line-buffer stream)))
      (do ((fa (read-fasta line-buffer)
               (read-fasta line-buffer))
           (counter 0 (1+ counter)))
          ((null fa) counter)))))

(defun benchmark-fasta2 (filename)
  (with-open-file (stream filename :direction :input
                   :element-type 'base-char
                   :external-format :ascii)
    (let ((line-buffer (make-line-buffer stream)))
      (do ((fa (read-fasta line-buffer)
               (read-fasta line-buffer))
           (counter 0 (1+ counter)))
          ((null fa) counter)))))



(defun split-fastq-file (filename chunk-size)
  (let ((file-pname (pathname filename)))
    (with-open-file (in file-pname :direction :input
                     :element-type '(unsigned-byte 8))
      (do* ((line-buffer (make-line-buffer in))
            (chunk-count 0 (1+ chunk-count))
            (chunk-pname (make-chunk-pname file-pname chunk-count)
                         (make-chunk-pname file-pname chunk-count)))
           ((not (more-lines-p line-buffer)))
        (write-n-fastq line-buffer chunk-size chunk-pname)))))



(defun write-n-fastq (line-buffer n chunk-pname)
  (with-open-file (out chunk-pname :direction :output
                   :if-exists :error
                   :element-type 'base-char)
    (do ((fq (read-fastq line-buffer #'write-fastq-sexp out)
             (read-fastq line-buffer #'write-fastq-sexp out))
         (count 1 (1+ count)))
        ((or (null fq)
             (= count n)) count))))

        


(defun benchmark2 (from to)
  (with-open-file (in from :direction :input
                   :element-type '(unsigned-byte 8))
    (with-open-file (out to :direction :output
                     :element-type 'base-char
                     :if-exists :supersede)
      (let ((line-buffer (make-line-buffer in)))
        (do ((fq (read-fastq line-buffer #'write-fastq-sexp out)
                 (read-fastq line-buffer #'write-fastq-sexp out))
             (counter 0 (1+ counter)))
            ((null fq) counter))))))

(defun benchmark3 (filename callback)
  (with-open-file (stream filename :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((line-buffer (make-line-buffer stream)))
      (loop
         for x = (read-fastq line-buffer callback :illumina)
         then (read-fastq line-buffer callback :illumina)
         while x
         collect x))))
    
(defun benchmark4 (filename)
  (with-open-file (stream filename :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((line-buffer (make-line-buffer stream)))
      (do ((fa (read-fasta line-buffer)
               (read-fasta line-buffer))
           (old-fa nil fa)
           (counter 0 (1+ counter)))
          ((null fa) (list counter
                           old-fa))))))

(defun benchmark5 (from to)
  (with-open-file (in from :direction :input
                   :element-type '(unsigned-byte 8))
    (with-open-file (out to :direction :output
                     :element-type 'base-char
                     :if-exists :supersede)
      (let ((line-buffer (make-line-buffer in)))
        (do ((fa (read-fasta line-buffer #'write-fasta-sexp out)
                 (read-fasta line-buffer #'write-fasta-sexp out))
             (counter 0 (1+ counter)))
            ((null fa) counter))))))
