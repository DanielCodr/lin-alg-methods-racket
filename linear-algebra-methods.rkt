;; The first three lines of this file were inserted by DrRacket. They record metadata
;; about the language level of this file in a form that our tools can easily process.
#reader(lib "htdp-intermediate-lambda-reader.ss" "lang")((modname linear-algebra-methods) (read-case-sensitive #t) (teachpacks ()) (htdp-settings #(#t constructor repeating-decimal #f #t none #f () #f)))
(define (gauss-jordan M)
  (top-down M empty))

(define (top-down M acc)
  (cond
    [(empty? M) (bottom-up acc empty)]
    [else
     (local
       [(define cols (columns M empty))
       (define nonzero-col (contains-nonzero cols 0))]
       (cond
         [(boolean? nonzero-col) (top-down empty acc)]
         [(zero? (first (list-ref cols nonzero-col))) (top-down (interchange M 0 (first-nonzero (list-ref cols nonzero-col) 0)) acc)]
         [(not (= 1 (first (list-ref cols nonzero-col))))
          (top-down (cons (map (λ (x) (/  x (leading-entry (first M)))) (first M)) (rest M)) acc)]
         [else (top-down-helper (first M) (first-nonzero (first M) 0) (leading-entry (first M)) (rest M) acc empty)]))]))

(define (top-down-helper ref-lst pivot-num pivot M acc reduced-acc)
  (cond
    [(empty? M) (top-down (reverse reduced-acc) (cons ref-lst acc))]
    [(zero? (list-ref (first M) pivot-num)) (top-down-helper ref-lst pivot-num pivot (rest M) acc (cons (first M) reduced-acc))]
    [else
     (top-down-helper ref-lst pivot-num pivot (rest M) acc
                      (cons (add (first M) (scale (* -1 (list-ref (first M) pivot-num)) ref-lst)) reduced-acc))]))

(define (top-down-helper2 ref-lst pivot-num pivot M acc reduced-acc)
  (cond
    [(empty? M) (reverse reduced-acc)]
    [(zero? (list-ref (first M) pivot-num)) (top-down-helper2 ref-lst pivot-num pivot (rest M) acc (cons (first M) reduced-acc))]
    [else
     (top-down-helper2 ref-lst pivot-num pivot (rest M) acc
                      (cons (add (first M) (scale (* -1 (list-ref (first M) pivot-num)) ref-lst)) reduced-acc))]))

(define (interchange M i1 i2)
  (interchange-helper M M i1 i2 0))
(define (interchange-helper M0 M i1 i2 count)
  (cond
    [(empty? M) empty]
    [(= count i1) (cons (list-ref M0 i2) (interchange-helper M0 (rest M) i1 i2 (add1 count)))]
    [(= count i2) (cons (list-ref M0 i1) (interchange-helper M0 (rest M) i1 i2 (add1 count)))]
    [else (cons (first M) (interchange-helper M0 (rest M) i1 i2 (add1 count)))]))

(define (add L1 L2)
  (add-helper L1 L2 empty))
(define (add-helper L1 L2 acc)
  (cond
    [(empty? L1) (reverse acc)]
    [else (add-helper (rest L1) (rest L2) (cons (+ (first L1) (first L2)) acc))]))
(define (scale k L)
  (map (λ (x) (* k x)) L))

(define (zero-row? r)
  (if (empty? r) true
      (if (zero? (first r)) (zero-row? (rest r))
          false)))

(define (column-helper M acc)
  (cond
    [(empty? M) (reverse acc)]
    [else (column-helper (rest M) (cons (first (first M)) acc))]))

(define (columns M acc)
  (cond
    [(empty? (first M)) (reverse acc)]
    [else (columns (map (λ (x) (rest x)) M) (cons (column-helper M empty) acc))]))

(define (contains-nonzero columns count)
  (cond
    [(empty? columns) false]
    [(zero-row? (first columns)) (contains-nonzero (rest columns) (+ count 1))]
    [else count]))

(define (first-nonzero col count)
  (cond
    [(zero? (first col)) (first-nonzero (rest col) (add1 count))]
    [else count]))


(define (leading-entry lst)
  (cond
    [(empty? lst) 'impossible]
    [(zero? (first lst)) (leading-entry (rest lst))]
    [else (first lst)]))
    
(define (bottom-up M acc)
  (cond
    [(empty? M) acc]
    [else (bottom-up (reduce (first M) (first-nonzero (first M) 0) (leading-entry (first M)) (rest M) empty) (cons (first M) acc))]))

(define (reduce ref-lst pivot-num pivot M acc)
  (cond
    [(empty? M) (reverse acc)]
    [(zero? (list-ref (first M) pivot-num)) (reduce ref-lst pivot-num pivot (rest M) (cons (first M) acc))]
    [else (reduce ref-lst pivot-num pivot (rest M) (cons (add (first M) (scale (* -1 (list-ref (first M) pivot-num)) ref-lst)) acc))]))

(define (identity-matrix n)
  (identity-matrix-helper n n empty))

(define (identity-matrix-helper sz n acc)
  (if (= n 0) acc
      (identity-matrix-helper sz (- n 1) (cons (place n sz 1) acc))))

(define (place n sz count)
  (if (= count n) (cons 1 (place n sz (+ count 1)))
      (if (> count sz) empty
          (cons 0 (place n sz (+ count 1))))))
(define (combine L1 L2)
  (if (empty? L1) empty
      (cons (append (first L1) (first L2)) (combine (rest L1) (rest L2)))))

(define (inverse A)
  (if (and (not (empty? A)) (= (length A) (length (first A))))
      (take (gauss-jordan (combine A (identity-matrix (length A)))) (length A))
      'impossible))

(define (take A m)
  (map (λ (x) (take-helper x m 0)) A))

(define (take-helper L m count)
  (if (empty? L) empty (if (< count m) (take-helper (rest L) m (+ count 1))
      (cons (first L) (take-helper (rest L) m (+ count 1))))))

(define (determinant A)
  (if (and (not (empty? A)) (= (length A) (length (first A))))
      (determinant-h A)
      'impossible))
(define (determinant-h  A)
  (if (= (length (first A)) 1) (first (first A))
      (determinant-helper (rest A) (first A) 1 0)))

(define (determinant-helper bot L count acc)
  (if (empty? L) acc
   (determinant-helper bot (rest L) (add1 count)
                       (+ acc (* (first L) (* (expt -1 (+ count 1)) (determinant-h (except bot count))))))))

(define (except A col)
  (map (λ (x) (except-helper x col 1)) A))

(define (except-helper L col count)
  (if (empty? L) empty
      (if (= col count) (except-helper (rest L) col (add1 count))
          (cons (first L) (except-helper (rest L) col (add1 count))))))

(define (matrix-mul A B)
  (map (λ (x) (map (λ (y) (dot x y)) (columns B empty))) A))

(define (matrix-vector-mul A v)
  (matrix-vector-mul-helper A v empty))

(define (matrix-vector-mul-helper A v acc)
  (if (empty? A) (reverse acc)
      (matrix-vector-mul-helper (rest A) v (cons (dot (first A) v) acc))))

(define (dot u v)
  (dot-helper u v 0))
(define (dot-helper u v acc)
  (if (empty? u) acc
      (dot-helper (rest u) (rest v) (+ acc (* (first u) (first v))))))
