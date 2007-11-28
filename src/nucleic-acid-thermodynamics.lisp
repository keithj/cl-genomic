;; Allawi and SantaLucia (1997) Biochemistry 36, 10581-10594
;;
;; 1M NaCl
;;
;; AA/TT, -7.9 ± 0.2, -22.2 ± 0.8
;; AT/TA, -7.2 ± 0.7, -20.4 ± 2.4
;; TA/AT, -7.2 ± 0.9, -21.3 ± 2.4
;; CA/GT, -8.5 ± 0.6, -22.7 ± 2.0
;; GT/CA, -8.4 ± 0.5, -22.4 ± 2.0
;; CT/GA, -7.8 ± 0.6, -21.0 ± 2.0
;; GA/CT, -8.2 ± 0.6, -22.2 ± 1.7
;; CG/GC, -10.6 ± 0.6, -27.2 ± 2.6
;; GC/CG, -9.8 ± 0.4, -24.4 ± 2.0
;; GG/CC, -8.0 ± 0.9, -19.9 ± 1.8
;; IGC, 0.1 ± 1.1, -2.8 ± 0.2
;; IAT, 2.3 ± 1.3, 4.1 ± 0.2

;; (probe-sequence target-sequence enthalpy entropy)
;; probe-sequence:  5'-3'
;; target-sequence: 5'-3'
;; enthalpy: (value error) delta-H in cal mol-1
;; entropy:  (value error) delta-S in cal mol-1 K-1

(defpackage na-therm
  (:use :common-lisp))

(in-package :na-therm)

(defconstant *molar-gas-constant* 1.987
  "Molar gas constant in cal mol-1.")
(defconstant *absolute-zero-celcius* -273.15
  "Zero Kelvin expressed in Celcius.")

(defvar *thermodynamic-parameters* 
  '((allawi-santalucia-1997
     ("AA" -7900 -22.2)
     ("AC" -8400 -22.4)
     ("AG" -7800 -21.0)
     ("AT" -7200 -20.4)
     ("CA" -8500 -22.7)
     ("CC" -8000 -19.9)
     ("CG" -10600 -27.2)
     ("CT" -7800 -21.0)
     ("GA" -8200 -22.2)
     ("GC" -9800 -24.4)
     ("GG" -8000 -19.9)
     ("GT" -8400 -22.4)
     ("TA" -7200 -21.3)
     ("TC" -8200 -22.2)
     ("TG" -8500 -22.7)
     ("TT" -7900 -22.2)
     ("G" 100 -2.8)
     ("C" 100 -2.8)
     ("A" 2300 4.1)
     ("T" 2300 4.1))))

(defvar *thermodynamic-parameter-tables* (make-hash-table))

(defun celcius-to-kelvin (celcius)
  (+ celcius (- *absolute-zero-celcius*)))

(defun kelvin-to-celcius (kelvin)
  (+ kelvin *absolute-zero-celcius*))

(defun get-thermodynamic-parameters (name)
  (unless (assoc name *thermodynamic-parameters*)
    (error "No such parameters as ~a~%" name))
  (multiple-value-bind (table cached)
      (gethash name *thermodynamic-parameter-tables*)
    (unless cached
      (setf table
            (setf (gethash name *thermodynamic-parameter-tables*)
                  (make-hash-table :test 'equal)))
      (dolist (datum (cdr (assoc name *thermodynamic-parameters*)))
        (setf (gethash (first datum) table) (rest datum))))
    table))

(defun order-2-with-termini (str)
  (let ((len (length str)))
    (if (zerop len)
        nil
      (append (cons (subseq str 0 1) (order-n str 2))
              (cons (subseq str (1- len) len) nil)))))

(defun order-n (str n)
  (let ((len (length str)))
    (if (< len n)
        nil
      (loop for i from 0 to (- len n)
            collect (subseq str i (+ i n))))))

(defun calculate-nn-parameters (str name)
  (let* ((tm-params (get-thermodynamic-parameters name))
         (param-pairs (mapcar #'(lambda (probe)
                                  (gethash probe tm-params))
                              (order-2-with-termini str))))
    (values (reduce #'+ param-pairs :key #'first)
            (reduce #'+ param-pairs :key #'second))))

(defun calculate-melting-temp (sequence seq-conc na-conc &key
                               (name 'allawi-santalucia-1997)
                               (f 4))
  (let ((phosphates (1- (length sequence))))
    (multiple-value-bind (enthalpy entropy)
        (calculate-nn-parameters sequence name)
      (kelvin-to-celcius (/ enthalpy
                            (+ (+ entropy (* 0.368 phosphates (log na-conc)))
                               (* *molar-gas-constant* (log (/ seq-conc f)))))))))

(defun read-sequences (filename)
  (with-open-file (stream filename)
    (do ((line (read-line stream nil)
               (read-line stream nil)))
        ((null line))
      (let ((first-tab (position #\tab line))
            (second-tab (position #\tab line :from-end t)))
        (format t "~a ~a ~a ~a~%"
                (subseq line 0 first-tab)
                (subseq line (1+ first-tab) second-tab)
                (subseq line second-tab)
                (calculate-melting-temp (subseq line (1+ first-tab) second-tab)
                                        0.000000001 0.05))))))
        