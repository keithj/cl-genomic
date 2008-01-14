
(in-package :bio-sequence)

(defparameter *fasta-line-width* 50)

(defmethod read-bio-sequence ((line-buffer byte-line-buffer) alphabet
                              (format (eql :fasta))
                              &optional (callback nil callback-supplied-p)
                              &rest callback-args)
  (let ((seq-header (find-line line-buffer #'fasta-header-p)))
    (if seq-header
        (multiple-value-bind (identity description)
            (parse-fasta-header (make-sb-string seq-header))
          (let ((seq-cache (make-array 0 :adjustable t :fill-pointer t)))
            (loop
               as line = (pull-line line-buffer)
               and cache-extend = (max 256 (floor (/ (length seq-cache) 2)))
               while line
               until (starts-with-byte-p line (char-code #\>))
               do (vector-push-extend line seq-cache cache-extend)
               finally (when line
                         (push-line line-buffer line)))
            (cond ((zerop (length seq-cache))
                   (error 'malformed-record-error :text
                          "incomplete Fasta record"))
                  (callback-supplied-p
                   (apply callback
                          (make-seq-sexp identity alphabet
                                         (concat-into-sb-string seq-cache)
                                         description) callback-args))
                  (t
                   (make-seq-sexp identity alphabet
                                  (concat-into-sb-string seq-cache)
                                  description)))))
      nil)))

(defun make-seq-fasta (sexp)
  "A callback which constructs a CLOS bio-sequence object from
sequence data that has been parsed into a standard SEXP."
  (let ((tag (car sexp))
        (content (cdr sexp)))      
    (unless (eql :bio-seq tag)
      (error "invalid :bio-seq sexp ~a" sexp))
    (make-seq :alphabet (assocdr :alphabet content)
              :ambiguity (assocdr :ambiguity content)
              :token-seq (assocdr :token-seq content)
              :name (assocdr :name content))))

(defun write-sexp-fasta (sexp &optional output-stream)
  "A callback which writes sequence data that has been parsed into a
standard SEXP to OUTPUT-STREAM in Fasta format."
  (if sexp
      (let* ((content (cdr sexp))
             (description (assocdr :description content)))
        (write-char #\> output-stream)
        (if (zerop (length description))
            (write-line (assocdr :identity content) output-stream)
          (progn
            (write-string (assocdr :identity content) output-stream)
            (write-char #\Space output-stream)
            (write-line description output-stream)))
        (write-wrapped-string (assocdr :token-seq content)
                              *fasta-line-width* output-stream))
    nil))

(defun fasta-header-p (bytes)
  (starts-with-byte-p bytes (char-code #\>)))

(defun parse-fasta-header (str)
  "Performs a basic parse of a Fasta header string STR by removing the
leading '>' character and splitting the line on the first space(s)
into identity and description. This function supports pathological
cases where the identity, description, or both are empty strings."
  (multiple-value-bind (split index)
      (split-sequence #\Space str :count 1 :remove-empty-subseqs t)
    (let ((str-len (length str))
          (str-elt-type (array-element-type str))
          (identity-str (car split)))
      (values 
       (if (> (length identity-str) 1)
           (adjust-array identity-str (- (length identity-str) 1)
                         :displaced-to identity-str
                         :displaced-index-offset 1)
         (make-string 0 :element-type str-elt-type))
       (if (< index str-len)
           (adjust-array str (- str-len index)
                         :displaced-to str
                         :displaced-index-offset index)
         (make-string 0 :element-type str-elt-type))))))
