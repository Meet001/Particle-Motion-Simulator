 #lang racket
(require "declarations.rkt")
;(require "drawing-routine.rkt")
;(require "testcases.rkt")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;e

(provide buildTree)
(provide calcForces)
(provide moveparticles)

(define (part-in? initialArea particle)
      (if (and (>= (vec-x (particle-posn particle)) (bbox-llx initialArea))
               (>= (vec-y (particle-posn particle)) (bbox-lly initialArea))
               (< (vec-x (particle-posn particle)) (bbox-rux initialArea))
               (< (vec-y (particle-posn particle)) (bbox-ruy initialArea)))#t
                                                                            #f))

(define (buildTree initialArea particles)
  (cond [(null? particles) (particle 0 (vec 0 0) (vec 0 0))]
        [(= (length particles) 1) (car particles)]
      [else (let*[(x1 (bbox-llx initialArea))
             (y1 (bbox-lly initialArea))
             (x2 (bbox-rux initialArea))
             (y2 (bbox-ruy initialArea))
             (x12 (/ (+ (bbox-llx initialArea) (bbox-rux initialArea)) 2))
             (y12 (/ (+ (bbox-lly initialArea) (bbox-ruy initialArea)) 2))
             (bbox1 (bbox x1 y1 x12 y12))
             (bbox2 (bbox x1 y12 x12 y2))
             (bbox3 (bbox x12 y1 x2 y12))
             (bbox4 (bbox x12 y12 x2 y2))]
        (gnode (total-mass particles) (centerof-mass particles)
               (list (buildTree bbox2 (part-inside bbox2 particles))
                     (buildTree bbox4 (part-inside bbox4 particles))
                     (buildTree bbox1 (part-inside bbox1 particles))
                     (buildTree bbox3 (part-inside bbox3 particles)))))]))
               
                    
  (define (part-inside initialArea particles)
  (define (help intialArea particles ans)
  (let* [(l (length particles))
        (x1 (bbox-llx initialArea))
        (y1 (bbox-lly initialArea))
        (x2 (bbox-rux initialArea))
        (y2 (bbox-ruy initialArea))]
    (cond [(null? particles) ans]
          [(part-in? initialArea (car particles))
           (help initialArea (cdr particles) (append ans (list (car particles))))]
          [else (help initialArea (cdr particles) ans)]
           )))
    (help initialArea particles '()))
          
        
(define (total-mass particles)
      (if (null? particles) 0
     (+ (particle-mass (car particles)) (total-mass (cdr particles)))))


(define (centerof-mass particles)
     (if (null? particles) 0
    (vec (/ (sum-mass-x particles) (total-mass particles)) (/ (sum-mass-y particles) (total-mass particles)))))


(define (sum-mass-x particles)
  (if (null? particles) 0
      (+ (* (vec-x (particle-posn (car particles))) (particle-mass (car particles))) (sum-mass-x (cdr particles)))))

(define (sum-mass-y particles)
  (if (null? particles) 0
      (+ (* (vec-y (particle-posn (car particles))) (particle-mass (car particles))) (sum-mass-y (cdr particles)))))




(define (calcForces initialArea tree particles)
   (let* [(x1 (bbox-llx initialArea))
        (y1 (bbox-lly initialArea))
        (x2 (bbox-rux initialArea))
        (y2 (bbox-ruy initialArea))
        (s (sqrt (/ (+ (sqr (- x2 x1)) (sqr (- y2 y1))) 2)))]
    (define (help tree particle side)
      (cond ((equal? tree particle) (vec 0 0))
            ((particle? tree) (force2 particle tree))
            ((= (length (gnode-subtrees tree)) 1) (force2 particle (car (gnode-subtrees tree))))
            ((> (/ (distance1 particle tree) side) 2)(force1 particle tree))
            (else (vec (+ (vec-x (help (list-ref (gnode-subtrees tree) 0) particle (/ side 2)))
                     (vec-x (help (list-ref (gnode-subtrees tree) 1) particle (/ side 2)))
                     (vec-x (help (list-ref (gnode-subtrees tree) 2) particle (/ side 2)))
                     (vec-x (help (list-ref (gnode-subtrees tree) 3) particle (/ side 2))))

                  (+ (vec-y (help (list-ref (gnode-subtrees tree) 0) particle (/ side 2)))
                     (vec-y (help (list-ref (gnode-subtrees tree) 1) particle (/ side 2)))
                     (vec-y (help (list-ref (gnode-subtrees tree) 2) particle (/ side 2)))
                     (vec-y (help (list-ref (gnode-subtrees tree) 3) particle (/ side 2))))))))
(if (null? particles) '()
    (append (list (help tree (car particles) s)) (calcForces initialArea tree (cdr particles))))))

(define (force2 particle1 particle2)
       (vec (* (/ (* 18 (particle-mass particle1) (particle-mass particle2))
                  (expt (distance2 particle1 particle2) 3))
               (- (vec-x (particle-posn particle2)) (vec-x (particle-posn particle1))))

            (* (/ (* 18 (particle-mass particle1) (particle-mass particle2))
                  (expt (distance2 particle1 particle2) 3))
               (- (vec-y (particle-posn particle2)) (vec-y (particle-posn particle1))))))



(define (force1 particle gnode)
       (vec (* (/ (* 18 (particle-mass particle) (gnode-mass gnode))
                  (expt (distance1 particle gnode) 3))
               (- (vec-x (gnode-posn gnode)) (vec-x (particle-posn particle))))

            (* (/ (* 18 (particle-mass particle) (gnode-mass gnode))
                  (expt (distance1 particle gnode) 3))
               (- (vec-y (gnode-posn gnode)) (vec-y (particle-posn particle))))))

(define (distance1 particle gnode)
 (sqrt (+ (sqr (- (vec-x (particle-posn particle)) (vec-x (gnode-posn gnode))))
  (sqr (- (vec-y (particle-posn particle)) (vec-y (gnode-posn gnode)))))))
    
(define (distance2 particle1 particle2)
 (sqrt (+ (sqr (- (vec-x (particle-posn particle1)) (vec-x (particle-posn particle2))))
  (sqr (- (vec-y (particle-posn particle1)) (vec-y (particle-posn particle2)))))))



(define (moveparticles particles forces)
  (define (help particle1 force)
      (particle (particle-mass particle1)
           (vec (+ (vec-x (particle-posn particle1)) (* (vec-x (particle-velocity particle1)) timeslice)
                (/ (* (vec-x force) timeslice timeslice) (* 2 (particle-mass particle1))))
                (+ (vec-y (particle-posn particle1)) (* (vec-y (particle-velocity particle1)) timeslice)
                (/ (* (vec-y force) timeslice timeslice) (* 2 (particle-mass particle1)))))
               
           (vec (+ (vec-x (particle-velocity particle1))
                   (/ (* (vec-x force) timeslice) (particle-mass particle1)))
                (+ (vec-y (particle-velocity particle1))
                   (/ (* (vec-y force) timeslice) (particle-mass particle1))))))
  (if (null? particles) '()
      (cons (help (car particles) (car forces)) (moveparticles (cdr particles) (cdr forces)))))






