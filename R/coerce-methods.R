setAs("GRanges", "VariantInfo",
      function(from)
      {
          class(from) <- "VariantInfo"
          mcols(from)$MAF <- as.numeric(rep(NA, length.out=length(from)))

          from
      })
