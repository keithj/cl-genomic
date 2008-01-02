
;; The car names the type of sexp, the cdr is an alist

(:bio-seq (:alphabet . :dna)
          (:name . "seq_name")
          (:description . "description line")
          (:ambiguity . nil)
          (:token-seq . "tgctagtcattcggatgcat"))

(:quality-seq (:bio-seq (:alphabet . :dna)
                        (:name . "seq_name")
                        (:description . "description line")
                        (:ambiguity . nil)
                        (:token-seq . "tgctagtcattcggatgcat"))
              (:quality . "<<<<<<<<<<<<<<<<<<<<"))

