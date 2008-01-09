
(in-package :bio-sequence)

(defparameter *seq-line-buffer-size* 512)

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

(defun make-chunk-pname (file-pname chunk-number)
  "Returns a new pathname for a file chunk based on a file pathname
FILE-PNAME and an integer CHUNK-NUMBER."
  (make-pathname
   :directory (pathname-directory file-pname)
   :name (concatenate 'string (pathname-name file-pname) "."
                      (princ-to-string chunk-number))
   :type (pathname-type file-pname)))


        
