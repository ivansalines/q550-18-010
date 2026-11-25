# Local validation — v18.53v anchor-alignment fix1

**Source repair only; no new physical IVPs and no full 204-node rerun.**

- Original v18.53v source manifest verification: PASS (14/14).
- Authenticated nested parents v18.53u (14 source, 8 output) and v18.53t (13 source, 7 output): PASS.
- Independent three-case matched-anchor replay: PASS for original/shared preformal metric at lambda 0, 1, 3.
- Original U archived arclength reconstructed exactly from its authenticated N10 arrays.
- Unaligned U/V arclength differences: -1.19278203571e-4, -2.02602351550e-4, -2.69815995098e-4 rad.
- Same-anchor U historical-algorithm vs V arclength differences: +4.27367798927e-8, -8.63846212340e-8, -8.79063399850e-9 rad.
- Original v18.53v ZIP packager regression test: PASS.
- Python syntax and bash syntax validation: PASS.
- Complete 36-row historical-anchor replay and 204-node output classification: **reserved for the user's server**, which has the original already-computed CSVs. The local tests cannot certify that classification.
