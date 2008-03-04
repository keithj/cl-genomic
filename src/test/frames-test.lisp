;;;
;;; Copyright (C) 2007-2008, Keith James. All rights reserved.
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

(in-package :cl-bio-test)

(fiveam:in-suite cl-bio-system:testsuite)


;;; Basic adding and removal of frames
(test add-frame/knowledgebase
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing")))
    (is (eql f (add-frame f kb)))
    (is-true (contains-frame-p f kb))
    ;; Adding a frame more than once is an error
    (signals knowledgebase-error
      (add-frame (make-instance 'frame :name "thing") kb))))

(test remove-frame/knowledgebase
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing")))
    (add-frame f kb)
    (is-true (contains-frame-p f kb))
    (is (eql f (remove-frame f kb)))
    (is-false (contains-frame-p f kb))
    (signals knowledgebase-error
      (remove-frame f kb))))


;;; Basic frame existence testing
(test contains-frame-p/knowledgebase
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing")))
    (add-frame f kb)
    (is-true (contains-frame-p f kb))))

(test find-frame/knowledgebase
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing")))
    (is-false (contains-frame-p "thing" kb))
    (add-frame f kb)
    (is (eql f (find-frame "thing" kb)))
    ;; Finding a non-existent frame is an error
    (signals knowledgebase-error
      (find-frame "no-such-frame" kb))))


;;; Basic adding and removal of slots
(test add-slot/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note")))
    (add-frame f kb)
    (is (eql f (add-slot f n)))
    (is (equalp (make-array 1 :initial-contents (list n)) (slots-of f)))
    ;; Adding a slot more than once is an error
    (signals knowledgebase-error
      (add-slot f n))))

(test remove-slot/knowledgebase
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note"))
        (s (make-instance 'single-valued-slot :name "score")))
    (add-frame f kb)
    ;; Test removal where slot is last in the slots vector
    (add-slot f n)
    (is-true (contains-slot-p f n))
    (is (eql f (remove-slot f n)))
    (is-false (contains-slot-p f n))
    ;; Test removal where slot is not last in the slots vector
    (add-slot f n)
    (add-slot f s)
    (is-true (contains-slot-p f n))
    (is-true (contains-slot-p f s))
    (is (eql f (remove-slot f s)))
    (is-true (contains-slot-p f n))
    (is-false (contains-slot-p f s))))

(test contains-slot-p/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note")))
    (add-frame f kb)
    (is-false (contains-slot-p f "note"))
    (is-false (contains-slot-p f n))
    (add-slot f n)
    (is-true (contains-slot-p f "note"))
    (is-true (contains-slot-p f n))))

(test find-slot/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note")))
    (add-frame f kb)
    (add-slot f n)
    (is (eql n (find-slot f "note")))))

(test value-of/frame/slot
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note")))
    (add-frame f kb)
    (add-slot f n)
    (is (null (value-of n)))
    (is (null (slot-value-of f "note")))
    (setf (slot-value-of f "note") "A note.")
    (is (string= "A note." (slot-value-of f "note")))))


;;; Adding a single-valued inverse slot


;;; Adding a set-valued inverse slot
(test value-of/frame/set-valued-inverse-slot
  (let ((kb (make-instance 'knowledgebase))
        (car (make-instance 'frame :name "car"))
        (wheel (make-instance 'frame :name "wheel"))
        (part-of-slot (make-instance 'part-of :name "part-of")))
    (add-frame car kb)
    (add-frame wheel kb)
    (add-slot wheel part-of-slot)
    (setf (slot-value-of wheel "part-of") (list car))
    (is-true (contains-slot-p car "has-part"))
    (is (equal (slot-value-of car "has-part") (list wheel)))))


;;; Adding a slot with a domain
(test add-slot/domain/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note" :domain 'frame)))
    (add-frame f kb)
    (is-false (contains-slot-p f n))
    (add-slot f n)
    (is-true (contains-slot-p f "note"))
    (signals knowledgebase-error
      (add-slot f (make-instance 'single-valued-slot :name "note"
                                 :domain 'simple-dna-sequence)))))

;;; Adding a slot with a range
(test add-slot/range/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'single-valued-slot :name "note" :range 'string)))
    (add-frame f kb)
    (is-false (contains-slot-p f "note"))
    (add-slot f n)
    (is (null (slot-value-of f "note")))
    (signals knowledgebase-error
      (setf (slot-value-of f "note") 999))
    (setf (slot-value-of f "note") "A note.")
    (is (string= "A note." (slot-value-of f "note")))))
