
(in-package :bio-sequence)

(defparameter *seq-line-buffer-size* 512)

(defun read-bio-sequence (line-input-stream &key alphabet ambiguity
                          format callback callback-args)
  "Reads a sequence record from LINE-INPUT-STREAM, optionally applying
function CALLBACK with additional CALLBACK-ARGS to the
result. Keywords are used to specify the expected alphabet (:dna,
:rna), ambiguity (:iupac, nil) and record format (:fasta, :fastq)."
  (unless alphabet
    (error "An alphabet must be supplied."))
  (unless format
    (error "A format must be supplied."))
  (when (and callback-args (not callback))
    (error "Callback arguments ~a were supplied without a callback."
           callback-args))
  (unless (member ambiguity '(:iupac :auto nil))
    (error "Invalid ambiguity (~a): expected one of ~a."
           ambiguity '(:iupac :auto nil)))
  (read-bio-sequence-alist line-input-stream format alphabet
                           ambiguity callback callback-args))

(defun make-seq-alist (identity alphabet ambiguity token-seq
                       &optional description)
  "Returns an alist, given a sequence IDENTITY, a vector of residue
tokens TOKEN-SEQ and a DESCRIPTION string."
  (pairlis '(:identity :alphabet :ambiguity :token-seq :description)
           (list identity alphabet
                 (cond ((and (eql :auto ambiguity)
                             (not (simplep token-seq alphabet)))
                        :iupac)
                       ((eql :auto ambiguity)
                        nil)
                       (t
                        ambiguity))
                 token-seq description)))

(defun make-quality-alist (identity alphabet ambiguity token-seq quality)
  "Returns an alist, given a sequence IDENTITY, a vector of residue
tokens TOKEN-SEQ and a QUALITY vector."
  (acons :quality quality
         (make-seq-alist identity alphabet ambiguity token-seq)))

(defun make-seq-from-alist (alist)
  "A callback which constructs a CLOS bio-sequence object from
sequence data that has been parsed into an ALIST."
  (make-seq  :alphabet (assocdr :alphabet alist)
             :ambiguity (assocdr :ambiguity alist)
             :identity (assocdr :identity alist)
             :token-seq (assocdr :token-seq alist)))

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

