
(in-package :bio-sequence)

(defmethod read-bio-sequence ((line-buffer byte-line-buffer) alphabet
                              (format (eql :fastq))
                              &optional (callback nil callback-supplied-p)
                              &rest callback-args)
  (let ((seq-header (loop
                       as line = (pull-line line-buffer)
                       while line
                       when (starts-with-byte-p line (char-code #\@))
                       return line)))
    (if seq-header
        (let ((seq (pull-line line-buffer))
              (quality-header (pull-line line-buffer))
              (quality (pull-line line-buffer)))
          (cond ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+))
                      callback-supplied-p)
                 (apply callback
                        (make-sexp-fastq seq-header
                                         seq alphabet quality)
                        callback-args))
                ((and seq quality-header quality
                      (starts-with-byte-p quality-header (char-code #\+)))
                 (make-sexp-fastq seq-header
                                  seq alphabet quality))
                (t
                 (error 'malformed-record-error :text
                        "incomplete Fastq record"))))
      nil)))

(defun make-quality-seq-fastq (sexp metric)
  "Callback which accepts a standard quality SEXP and creates a new
dna-quality-sequence with quality METRIC."
  (if sexp
      (let ((content (cdr sexp)))
        (make-quality-seq :alphabet (assocdr :alphabet content)
                          :ambiguity (assocdr :ambiguity content)
                          :token-seq (assocdr :token-seq content)
                          :quality (assocdr :quality content)
                          :name (assocdr :name content)
                          :metric metric))
    nil))

(defun write-sexp-fastq (sexp &optional output-stream)
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
                                #'write-sexp-fastq out)
             (read-bio-sequence line-buffer :dna :fastq
                                #'write-sexp-fastq out))
         (count 1 (1+ count)))
        ((or (null fq)
             (= count n)) count))))

(defun make-sexp-fastq (seq-header seq alphabet quality)
  "Creates a canonical quality sexp from the raw byte arrays of the
header SEQ-HEADER (which must include the '@' character), the sequence
SEQ and base quality QUALITY."
  (make-quality-sexp (make-seq-sexp (make-sb-string seq-header 1)
                                    (make-sb-string seq))
                     alphabet
                     (make-sb-string quality)))