
(require :sb-cover)

(declaim (optimize sb-cover:store-coverage-data))

(asdf:operate 'asdf:load-op :cl-genomic)
(asdf:operate 'asdf:test-op :cl-genomic)

(sb-cover:report "./test-coverage-report/")

(declaim (optimize (sb-cover:store-coverage-data 0)))

(asdf:operate 'asdf:compile-op :cl-genomic)

(quit)
