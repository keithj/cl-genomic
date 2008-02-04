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
    (signals knowledgebase-error
      (find-frame "no-such-frame" kb))))

;;; Basic adding and removal of slots
(test add-slot/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'slot :name "note")))
    (signals knowledgebase-error
      (add-slot f n kb))
    (add-frame f kb)
    (is (eql f (add-slot f n kb)))
    (is (equalp (make-array 1 :initial-contents (list n)) (slots-of f)))
    (signals knowledgebase-error
      (add-slot f n kb))))

(test remove-slot/knowledgebase
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'slot :name "note"))
        (s (make-instance 'slot :name "score")))
    (add-frame f kb)
    ;; Test removal where slot is last in the slots vector
    (add-slot f n kb)
    (is-true (contains-slot-p f n kb))
    (is (eql f (remove-slot f n kb)))
    (is-false (contains-slot-p f n kb))
    ;; Test removal where slot id not last in the slots vector
    (add-slot f n kb)
    (add-slot f s kb)
    (is-true (contains-slot-p f n kb))
    (is-true (contains-slot-p f s kb))
    (is (eql f (remove-slot f s kb)))
    (is-true (contains-slot-p f n kb))
    (is-false (contains-slot-p f s kb))))

(test contains-slot-p/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'slot :name "note")))
    (add-frame f kb)
    (is-false (contains-slot-p f "note" kb))
    (is-false (contains-slot-p f n kb))
    (add-slot f n kb)
    (is-true (contains-slot-p f "note" kb))
    (is-true (contains-slot-p f n kb))))

(test find-slot/frame
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'slot :name "note")))
    (add-frame f kb)
    (add-slot f n kb)
    (is (eql n (find-slot f "note" kb)))))

(test value-of/frame/slot
  (let ((kb (make-instance 'knowledgebase))
        (f (make-instance 'frame :name "thing"))
        (n (make-instance 'slot :name "note")))
    (add-frame f kb)
    (add-slot f n kb)
    (is (null (value-of n)))
    (is (null (slot-value-of f n kb)))
    (is (null (slot-value-of "thing" "note" kb)))
    (setf (slot-value-of f n kb) "A note.")
    (is (string= "A note." (slot-value-of f n kb)))))

(test value-of/frame/inverse-slot
  (let ((kb (make-instance 'knowledgebase))
        (car (make-instance 'frame :name "car"))
        (wheel (make-instance 'frame :name "wheel"))
        (wheel-slot (make-instance 'part-of :name "car-part")))
    (add-frame car kb)
    (add-frame wheel kb)
    (add-slot car wheel-slot kb)
    (setf (slot-value-of car wheel-slot kb) wheel)))