
(require :sb-cover)

(declaim (optimize sb-cover:store-coverage-data))

(asdf:oos 'asdf:load-op :cl-genomic)

(asdf:oos 'asdf:test-op :cl-gp-utilities)
(asdf:oos 'asdf:test-op :cl-io-utilities)
(asdf:oos 'asdf:test-op :cl-genomic)

(sb-cover:report "./test-coverage-report/")

(declaim (optimize (sb-cover:store-coverage-data 0)))

(asdf:oos 'asdf:compile-op :cl-genomic :force t)

(quit)
