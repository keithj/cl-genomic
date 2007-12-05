
(in-package :bio-sequence)

(enable-lazy-init (find-class 'simple-dna-quality-sequence))

(defun as-string (bytes)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type (simple-array (unsigned-byte 8)) bytes))
  (let ((string (make-array (length bytes)
                            :element-type 'base-char)))
    (declare (type simple-base-string string))
    (copy-array bytes 0 (1- (length bytes))
                string 0 #'code-char)
    string))

(defmethod read-fastq ((obj line-buffer))
;;   (loop as line = (pull-line obj)
;;      while line
;;      do (when (and (not (zerop (length line)))
;;                    (= (char-code #\@) (aref line 0)))
;;           (push-line obj line)))
  (declare (optimize (speed 3) (safety 1)))
  (let ((seq-header (pull-line obj))
        (seq (pull-line obj))
        (qual-header (pull-line obj))
        (quality (pull-line obj)))
    (declare (type (simple-array (unsigned-byte 8)) seq))
    (cond ((and seq-header (or (position (char-code #\N) seq)
                               (position (char-code #\n) seq)))
           (make-instance 'lazy-init-proxy
                          :proxied-class 'iupac-dna-quality-sequence
                          :initargs
                          (list :token-seq (encode-iupac-seq
                                            (as-string seq)
                                            #'encode-dna-4bit)
                                :quality (as-string quality))))
          (seq-header
           (make-instance 'lazy-init-proxy
                          :proxied-class 'simple-dna-quality-sequence
                          :initargs
                          (list :token-seq (encode-simple-seq
                                            (as-string seq)
                                            #'encode-dna-2bit)
                                :quality (as-string quality))))
          (t
           nil))))


(defun benchmark (filename)
  (with-open-file (stream filename :direction :input
                   :element-type '(unsigned-byte 8))
    (let ((line-buffer (make-line-buffer stream)))
      (do ((fq (read-fastq line-buffer)
               (read-fastq line-buffer))
           (old-fq nil fq)
           (counter 0 (1+ counter)))
          ((null fq) (list counter
                           old-fq
                           (token-seq-of old-fq)
                           old-fq))))))
