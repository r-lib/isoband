# basic functions

    Code
      isolines_grob(l, labels = c("a", "b", "c"))
    Condition
      Error in `isolines_grob()`:
      ! Number of labels must match the number of breaks.

---

    Code
      isolines_grob(l, margin = 1:4)
    Condition
      Error in `isolines_grob()`:
      ! The `margin` parameter must be a unit object of length four.

---

    Code
      isolines_grob(l, margin = grid::unit(1:3, "pt"))
    Condition
      Error in `isolines_grob()`:
      ! The `margin` parameter must be a unit object of length four.

