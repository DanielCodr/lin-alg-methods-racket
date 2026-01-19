This racket module provides core linear algebra functionality, including Gaussian elimination, Matrix inverses, determinants, and matrix and vector operations.

We list some example usage below.

(define A '((1 2 3) (0 1 4) (5 6 0)))

;; Row-reduce A
(gauss-jordan A)

;; Compute inverse
(inverse A)

;; Compute determinant
(determinant A)

;; Multiply matrices
(matrix-mul A (inverse A)) ;; should return identity matrix

;; Multiply matrix by vector
(matrix-vector-mul A '(1 0 0)) ;; returns first column of A

;; Dot product
(dot '(1 2 3) '(4 5 6)) ;; returns 32
