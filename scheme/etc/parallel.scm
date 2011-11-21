(define (make n)
  (letrec (((iter i acc) (if (>= i n)
                             acc
                             (iter (+ i 1) (cons i acc)))))
    (iter 0 '())))

(define (sum l) (apply + l))

(define (even? n) (= (* (/ n 2) 2) n))

(definerec (map f l)
  (if (null? l)
      '()
      (cons (f (car l)) (map f (cdr l)))))

(definerec (collatz n)
  (cond ((<= n 4) '())
        ((even? n) (cons n (collatz (/ n 2))))
        (else (cons n (collatz (+ (* 3 n) 1))))))

(let ((x (parallel-map collatz (make 10000))))
  (display "done!"))