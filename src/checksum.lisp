;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
;;;
;;; This file is part of cl-genomic.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :bio-sequence)

(defmethod seguid ((seq encoded-dna-sequence) &key (start 0) end)
  (base64:usb8-array-to-base64-string
   (digest-encoded-sequence seq #'decode-dna-4bit :sha1 start end)))

(defmethod seguid ((seq encoded-rna-sequence) &key (start 0) end)
  (base64:usb8-array-to-base64-string
   (digest-encoded-sequence seq #'decode-rna-4bit :sha1 start end)))

(defmethod seguid ((seq encoded-aa-sequence) &key (start 0) end)
  (base64:usb8-array-to-base64-string
   (digest-encoded-sequence seq #'decode-aa-7bit :sha1 start end)))

(defmethod seguid ((seq mapped-dna-sequence) &key (start 0) end)
  (base64:usb8-array-to-base64-string
   (digest-mapped-sequence seq :sha1 start end)))


(defmethod md5 ((seq encoded-dna-sequence) &key (start 0) end)
  (ironclad:byte-array-to-hex-string 
   (digest-encoded-sequence seq #'decode-dna-4bit :md5 start end)))

(defmethod md5 ((seq encoded-rna-sequence) &key (start 0) end)
  (ironclad:byte-array-to-hex-string
   (digest-encoded-sequence seq #'decode-rna-4bit :md5 start end)))

(defmethod md5 ((seq encoded-aa-sequence) &key (start 0) end)
  (ironclad:byte-array-to-hex-string
   (digest-encoded-sequence seq #'decode-aa-7bit :md5 start end)))

(defmethod md5 ((seq mapped-dna-sequence) &key (start 0) end)
  (ironclad:byte-array-to-hex-string
   (digest-mapped-sequence seq :md5 start end)))

(defun digest-encoded-sequence (seq decoder digest start &optional end)
  (with-accessors ((vector vector-of))
      seq
    (let ((end (or end (length vector)))
          (stream (crypto:make-digesting-stream digest)))
      (loop
         for i of-type fixnum from start below end
         do (stream-write-byte stream (char-code
                                       (char-upcase
                                        (funcall decoder (aref vector i))))))
       (ironclad:produce-digest stream))))

(defun digest-mapped-sequence (seq digest start &optional end)
  (declare (optimize (speed 3)))
  (with-slots ((length dxn:length) (area dxn:mmap-area))
      seq
    (let ((end (or end length))
          (ptr (dxn:mmap-area-ptr area)))
      (if (dxn:mmap-area-live-p area) ; Only access residues if mapped
          (let ((stream (crypto:make-digesting-stream digest)))
            (loop
               for i of-type fixnum from start below end
               do (stream-write-byte
                   stream (char-code
                           (char-upcase
                            (code-char (cffi:mem-aref ptr :char i))))))
             (ironclad:produce-digest stream))))))
