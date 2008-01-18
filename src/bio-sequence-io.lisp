
(in-package :bio-sequence)

(defparameter *seq-line-buffer-size* 512)

(defun read-bio-sequence (line-buffer &key alphabet ambiguity format
                          callback callback-args)
  "Reads a sequence record from LINE-BUFFER, optionally applying
function CALLBACK with additional CALLBACK-ARGS to the
result. Keywords are used to specify the expected alphabet (:dna,
:rna), ambiguity (:iupac, nil) and record format (:fasta, :fastq)."
  (unless alphabet
    (error "an alphabet must be supplied"))
  (unless format
    (error "a format must be supplied"))
  (when (and callback-args
             (not callback))
    (error "callback-args (~a) were supplied without a callback"
           callback-args))
  (unless (member ambiguity '(:iupac :auto nil))
    (error "invalid ambiguity (~a), expected one of ~a"
           ambiguity '(:iupac :auto nil)))
  (read-bio-sequence-alist line-buffer format alphabet ambiguity
                           callback callback-args))

(defun make-seq-alist (name alphabet ambiguity token-seq
                       &optional description)
  "Returns an alist, given a sequence NAME, a vector of residue tokens
TOKEN-SEQ and a DESCRIPTION string."
  (pairlis '(:name :alphabet :ambiguity :token-seq :description)
           (list name alphabet
                 (cond ((and (eql :auto ambiguity)
                             (not (simplep token-seq alphabet)))
                        :iupac)
                       ((eql :auto ambiguity)
                        nil)
                       (t
                        ambiguity))
                 token-seq description)))

(defun make-quality-alist (name alphabet ambiguity token-seq quality)
  "Returns an alist, given a sequence NAME, a vector of residue tokens
TOKEN-SEQ and a QUALITY vector."
  (acons :quality quality
         (make-seq-alist name alphabet ambiguity token-seq)))

(defun make-seq-from-alist (alist)
  "A callback which constructs a CLOS bio-sequence object from
sequence data that has been parsed into an ALIST."
  (make-seq :alphabet (assocdr :alphabet alist)
            :ambiguity (assocdr :ambiguity alist)
            :token-seq (assocdr :token-seq alist)
            :name (assocdr :name alist)))

(defun make-removing-callback (callback predicate)
  (lambda (sexp &rest callback-args)
    (if (funcall predicate sexp)
        nil
      (apply callback sexp callback-args))))

(defun make-chunk-pname (file-pname chunk-number)
  "Returns a new pathname for a file chunk based on a file pathname
FILE-PNAME and an integer CHUNK-NUMBER."
  (make-pathname
   :directory (pathname-directory file-pname)
   :name (concatenate 'string (pathname-name file-pname) "."
                      (princ-to-string chunk-number))
   :type (pathname-type file-pname)))

