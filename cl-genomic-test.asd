;;;
;;; Copyright (C) 2007-2010 Keith James. All rights reserved.
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

(defsystem cl-genomic-test
    :depends-on (:cl-genomic
                 (:version :deoxybyte-utilities "0.6.0")
                 (:version :deoxybyte-io "0.6.0")
                 (:version :lift "1.7.0"))
    :components ((:module :cl-genomic-test
                          :pathname "test/"
                          :components
                          ((:file "package")
                           (:file "cl-genomic-test"
                                  :depends-on ("package"))
                           (:file "bio-sequence-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"))
                           (:file "bio-sequence-io-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"
                                               "bio-sequence-test"))
                           (:file "bio-sequence-encoding-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"
                                               "bio-sequence-test"))
                           (:file "bio-sequence-translation-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"))
                           (:file "bio-sequence-interval-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"))
                           (:file "bio-sequence-alignment-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"))
                           (:file "fjoin-test"
                                  :depends-on ("package"
                                               "cl-genomic-test"))))))
