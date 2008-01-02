
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
                        (make-fastq-sexp seq-header seq quality)
                        callback-args))
                ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+)))
                 (make-fastq-sexp seq-header seq quality))
                (t
                 (error 'malformed-record-error
                        :text "incomplete Fastq record"))))
      nil)))

(defun make-fastq-sexp (seq-header seq quality)
  "Creates a canonical Fastq sexp from the raw byte arrays of the
header (including the '@' character) SEQ-HEADER, the sequence SEQ and
base quality QUALITY."
  (make-quality-sexp (make-dna-sexp (make-sb-string seq-header 1)
                                    (make-sb-string seq))
                     (make-sb-string quality)))

(defun write-fastq-sexp (sexp &optional output-stream)
  "Callback which accepts a canonical Fastq SEXP and writes it to
OUTPUT-STREAM as a Fastq format record. OUTPUT-STREAM defaults to
*standard-output*."
  (if sexp
      (let ((seq (cdadr sexp))
            (qual (cddr sexp)))
        (write-line (assocdr :name seq) output-stream)
        (write-line (assocdr :token-seq seq) output-stream)
        (write-line "+" output-stream)
        (write-line (assocdr :quality qual) output-stream)
        t)
    nil))

(defun make-fastq-seq (sexp metric)
  "Callback which accepts a standard Fastq SEXP and creates a new
dna-quality-sequence with quality METRIC."
  (if sexp
      (let ((content (cdr sexp)))
        (make-dna-quality-seq (assocdr :token-seq content)
                              (assocdr :quality content)
                              :ambiguity (assocdr :ambiguity content)
                              :metric metric))
    nil))
