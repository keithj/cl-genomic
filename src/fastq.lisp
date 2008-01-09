
(in-package :bio-sequence)

(defmethod read-fastq ((obj byte-line-buffer) &optional
                       (callback nil callback-supplied-p)
                       &rest callback-args)
  (let ((seq-header (loop
                       as line = (pull-line obj)
                       while line
                       when (starts-with-byte-p line (char-code #\@))
                       return line)))
    (if seq-header
        (let ((seq (pull-line obj))
              (quality-header (pull-line obj))
              (quality (pull-line obj)))
          (cond ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+))
                      callback-supplied-p)
                 (apply callback
                        (make-quality-sexp-fastq seq-header seq quality)
                        callback-args))
                ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+)))
                 (make-quality-sexp-fastq seq-header seq quality))
                (t
                 (error 'malformed-record-error :text
                        "incomplete Fastq record"))))
      nil)))

(defun make-quality-sexp-fastq (seq-header seq quality)
  "Creates a canonical quality sexp from the raw byte arrays of the
header SEQ-HEADER (which must include the '@' character), the sequence
SEQ and base quality QUALITY."
  (make-quality-sexp (make-dna-sexp (make-sb-string seq-header 1)
                                    (make-sb-string seq))
                     (make-sb-string quality)))

(defun write-qual-sexp-fastq (sexp &optional output-stream)
  "Callback which accepts a canonical quality SEXP and writes it to
OUTPUT-STREAM as a Fastq format record. OUTPUT-STREAM defaults to
*standard-output*."
  (if sexp
      (let ((seq (cdadr sexp))
            (qual (cddr sexp)))
        (write-char #\@ output-stream)
        (write-line (assocdr :name seq) output-stream)
        (write-line (assocdr :token-seq seq) output-stream)
        (write-line "+" output-stream)
        (write-line (assocdr :quality qual) output-stream)
        t)
    nil))

(defun make-qual-seq-from-sexp (sexp metric)
  "Callback which accepts a standard quality SEXP and creates a new
dna-quality-sequence with quality METRIC."
  (if sexp
      (let ((content (cdr sexp)))
        (make-dna-quality-seq (assocdr :token-seq content)
                              (assocdr :quality content)
                              :ambiguity (assocdr :ambiguity content)
                              :metric metric))
    nil))

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
    (do ((fq (read-fastq line-buffer #'write-qual-sexp-fastq out)
             (read-fastq line-buffer #'write-qual-sexp-fastq out))
         (count 1 (1+ count)))
        ((or (null fq)
             (= count n)) count))))
