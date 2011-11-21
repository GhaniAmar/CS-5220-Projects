; Test of mutual recursion with even?
(define mr-even?
  (letrec (((even x)
            (if (= x 0) #t (odd (- x 1))))
           ((odd x)
            (if (= x 0) #f (even (- x 1)))))
          even))

; Ackermann's function
(definerec (ackermann x y)
  (cond ((= y 0) 0)
        ((= x 0) (* 2 y))
        ((= y 1) 2)
        (else (ackermann (- x 1)
                         (ackermann x (- y 1))))))

; Factorial using apply keyword
(define (factorial n)
  (letrec (((niter i)
            (if (> i n)
                '()
                (cons i (niter (+ i 1))))))
    (apply * (niter 1))))

(define (even? n)
  (= (* (/ n 2) 2) n))

; Logarithmic time exponentiation
(define (expt base exponent)
  (letrec (((iter-exp a b n)
            (cond ((= n 0) a)
                  ((even? n) (iter-exp a (* b b) (/ n 2)))
                  (else (iter-exp (* a b) (* b b) (/ (- n 1) 2))))))
    (iter-exp 1 base exponent)))

(define (square x) (* x x))

; Logarithmic time nth Fibonacci number
(define (fib n)
  (letrec (((fib-iter a b p q count)
            (cond ((= count 0) b)
                  ((even? count)
                   (fib-iter a
                             b
                             (+ (square p) (square q))
                             (+ (* 2 p q) (square q))
                             (/ count 2)))
                  (else (fib-iter (+ (* b q) (* a q) (* a p)) (+ (* b p) (* a q))
                                  p q
                                  (- count 1))))))
    (fib-iter 1 0 0 1 n)))

; Church numerals
(define zero (lambda (f) (lambda (x) x)))
(define (succ n)
  (lambda (f) (lambda (x) (f ((n f) x)))))
(define (add x y) ((x succ) y))
(define (comp-add n m)
  (lambda (f) (lambda (x) ((n f) ((m f) x)))))
(define (to-int n)
  ((n (lambda (x) (+ 1 x))) 0))
(definerec (from-int n)
  (if (= n 0)
      zero
      (succ (from-int (- n 1)))))
