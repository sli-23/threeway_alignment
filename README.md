# Three-Way Sequence Alignment
Design a **dynamic programming algorithm** to align 3 sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub> such that **the sum of pairwise alignment scores** is maximized, given the **BLOSUM62**<sup>1</sup> substitution model and **non-affine gap penalty**.

### Input
The `threeway_align` function takes as input the following parameters:
* Three Python [strings](https://docs.python.org/3/library/stdtypes.html#textseq) `s1`, `s2`, and `s3`, representing the three amino acid sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub>
* A real-valued number (usually is negative) `gap`, which is the gap penalty.

### Output
The `threeway_align` function returns three strings `aligned_s1`, `aligned_s2`, and `aligned_s3` with the following properties:
* All three have equal length
    * `len(aligned_s1) == len(aligned_s2) == len(aligned_s3)`
* Removing all gap characters from `aligned_s1`, `aligned_s2`, and `aligned_s3` yields `s1`, `s2`, and `s3`, respectively
    * `aligned_s1.replace('-','') == s1 and aligned_s2.replace('-','') == s2 and aligned_s3.replace('-','') == s3`
