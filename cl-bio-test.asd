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

(in-package #:cl-bio-system)

(defsystem cl-bio-test
    :depends-on (:cl-bio :fiveam)
    :components ((:module :cl-bio-test
                          :pathname "src/test/"
                          :components
                          ((:file "package")
                           (:file "bio-sequence-test"
                                  :depends-on ("package"))
                           (:file "bio-sequence-io-test"
                                  :depends-on ("package"
                                               "bio-sequence-test"))
                           (:file "frames-test"
                                  :depends-on ("package"))))))
