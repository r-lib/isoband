# Inconsistent numbers of isoband levels cause an error

    Code
      isobands(1:4, 1:4, m, c(0.5, 1.5, 2.5), c(0.5, 1.5))
    Condition
      Error in `isobands()`:
      ! Vectors specifying isoband levels must be of equal length or of length 1

---

    Code
      isobands(1:4, 1:4, m, c(0.5, 1.5), c(0.5, 1.5, 2.5))
    Condition
      Error in `isobands()`:
      ! Vectors specifying isoband levels must be of equal length or of length 1

