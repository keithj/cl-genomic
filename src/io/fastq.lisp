
(in-package :bio-sequence)

(defmethod read-bio-sequence-alist ((line-buffer byte-line-buffer) (format (eql :fastq))
                                    alphabet ambiguity
                                    &optional (callback nil callback-supplied-p)
                                    callback-args)
  (let ((seq-header (find-line line-buffer #'fastq-header-p)))
    (if seq-header
        (let ((seq (pull-line line-buffer))
              (quality-header (pull-line line-buffer))
              (quality (pull-line line-buffer)))
          (cond ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+))
                      callback-supplied-p)
                 (apply callback (make-quality-alist (make-sb-string seq-header 1)
                                                     alphabet ambiguity
                                                     (make-sb-string seq)
                                                     quality)
                        callback-args))
                ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+)))
                 (make-quality-alist (make-sb-string seq-header 1)
                                     alphabet ambiguity
                                     (make-sb-string seq)
                                     quality))
                (t
                 (error 'malformed-record-error :text
                        "incomplete Fastq record"))))
      nil)))

(defun make-quality-seq-fastq (alist metric)
  "Callback which accepts a ALIST and creates a new
dna-quality-sequence with quality METRIC."
  (make-quality-seq :alphabet (assocdr :alphabet alist)
                    :ambiguity (assocdr :ambiguity alist)
                    :token-seq (assocdr :token-seq alist)
                    :quality (assocdr :quality alist)
                    :name (assocdr :name alist)
                    :metric metric))

(defun write-alist-fastq (alist &optional output-stream)
  "Callback which accepts an ALIST and writes it to OUTPUT-STREAM as a
Fastq format record. OUTPUT-STREAM defaults to *standard-output*."
  (write-char #\@ output-stream)
  (write-line (assocdr :name alist) output-stream)
  (write-line (assocdr :token-seq alist) output-stream)
  (write-line "+" output-stream)
  (write-line (assocdr :quality alist) output-stream)
  t)

(defun split-fastq-file (filename chunk-size)
  "Splits Fastq file FILENAME into automatically named chunks, each,
except the last file, containing up to CHUNK-SIZE records."
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
  "Reads up to N Fastq records from LINE-BUFFER and writes them into a
new file of pathname CHUNK-PNAME."
  (with-open-file (out chunk-pname :direction :output
                   :if-exists :error
                   :element-type 'base-char)
    (do ((fq (read-bio-sequence line-buffer :dna :fastq
                                #'write-alist-fastq out)
             (read-bio-sequence line-buffer :dna :fastq
                                #'write-alist-fastq out))
         (count 1 (1+ count)))
        ((or (null fq)
             (= count n)) count))))

(defun fastq-header-p (bytes)
  (starts-with-byte-p bytes (char-code #\@)))
