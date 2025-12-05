# empty or incorrect input

    Code
      clip_lines_impl(numeric(0), numeric(1), integer(0), 3, 2, 0.1, 0.1, 0)
    Condition
      Error in `clip_lines_impl()`:
      ! Number of x and y coordinates must match.

---

    Code
      clip_lines_impl(numeric(0), numeric(0), integer(1), 3, 2, 0.1, 0.1, 0)
    Condition
      Error in `clip_lines_impl()`:
      ! Number of x coordinates and id values must match.

---

    Code
      clip_lines_impl(numeric(1), numeric(0), integer(0), 3, 2, 0.1, 0.1, 0)
    Condition
      Error in `clip_lines_impl()`:
      ! Number of x and y coordinates must match.

