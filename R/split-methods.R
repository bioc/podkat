setMethod("split", signature=signature(x="GRanges", f="GRangesList"),
          function(x, f)
          {
              out <- GRangesList(lapply(f, function(y) intersect(x, y)))
              mcols(out) <- mcols(f)
              out[which(sapply(out, length) > 0)]
          })
