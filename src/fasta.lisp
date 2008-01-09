
(in-package :bio-sequence)

(defparameter *fasta-line-width* 50)

(defmethod read-fasta ((obj byte-line-buffer) &optional
                       (callback nil callback-supplied-p)
                       &rest callback-args)
  (let ((seq-header (loop as line = (pull-line obj)
                       while line
                       when (starts-with-byte-p line (char-code #\>))
                       return line)))
    (if seq-header
        (multiple-value-bind (name description)
            (parse-fasta-header (make-sb-string seq-header))
          (let ((seq-cache (make-array 0 :adjustable t :fill-pointer t)))
            (loop
               as line = (pull-line obj)
               and cache-extend = (max 256 (floor (/ (length seq-cache) 2)))
               while line
               until (starts-with-byte-p line (char-code #\>))
               do (vector-push-extend line seq-cache cache-extend)
               finally (when line
                         (push-line obj line)))
            (cond ((zerop (length seq-cache))
                   (error 'malformed-record-error :text
                          "incomplete Fasta record"))
                  (callback-supplied-p
                   (apply callback
                          (make-dna-sexp name
                                         (concat-into-sb-string seq-cache)
                                         description) callback-args))
                  (t
                   (make-dna-sexp name (concat-into-sb-string seq-cache)
                                  description)))))

      nil)))

(defmethod read-fasta ((obj line-buffer) &optional
                       (callback nil callback-supplied-p)
                       &rest callback-args)
  (let ((seq-header (loop as line = (pull-line obj)
                       while line
                       when (starts-with-char-p line #\>)
                       return line)))
    (if seq-header
        (multiple-value-bind (name description)
            (parse-fasta-header seq-header)
          (let ((seq-cache (make-array 0 :adjustable t :fill-pointer t)))
            (loop
               as line = (pull-line obj)
               and cache-extend = (max 256 (floor (/ (length seq-cache) 2)))
               while line
               until (starts-with-char-p line #\>)
               do (vector-push-extend line seq-cache cache-extend)
               finally (when line
                         (push-line obj line)))
            (cond ((zerop (length seq-cache))
                   (error 'malformed-record-error :text
                          "incomplete Fasta record"))
                  (callback-supplied-p
                   (apply callback
                          (make-dna-sexp name
                                         (concat-strings seq-cache)
                                         description) callback-args))
                  (t
                   (make-dna-sexp name (concat-strings seq-cache)
                                  description)))))

      nil)))

(defun parse-fasta-header (str)
  "Performs a basic parse of a Fasta header string STR by removing the
leading '>' character and splitting the line on the first space(s)
into name and description. This function supports pathological cases
where the name, description, or both are empty strings."
  (multiple-value-bind (split index)
      (split-sequence #\Space str :count 1 :remove-empty-subseqs t)
    (let ((str-len (length str))
          (str-elt-type (array-element-type str))
          (name-str (car split)))
      (values 
       (if (> (length name-str) 1)
           (adjust-array name-str (1- str-len)
                         :displaced-to name-str
                         :displaced-index-offset 1)
         (make-string 0 :element-type str-elt-type))
       (if (< index str-len)
           (adjust-array str (- str-len index)
                         :displaced-to str
                         :displaced-index-offset index)
         (make-string 0 :element-type str-elt-type))))))

(defun write-dna-sexp-fasta (sexp &optional output-stream)
  (if sexp
      (let* ((content (cdr sexp))
             (description (assocdr :description content)))
        (write-char #\> output-stream)
        (if (zerop (length description))
            (write-line (assocdr :name content) output-stream)
          (progn
            (write-string (assocdr :name content) output-stream)
            (write-char #\Space output-stream)
            (write-line description output-stream)))
        (write-wrapped-string (assocdr :token-seq content)
                              *fasta-line-width* output-stream))
    nil))
