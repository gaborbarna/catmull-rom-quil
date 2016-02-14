(ns catmull-rom-quil.core
  (:refer-clojure :exclude [* - + == / < <= > >= not= = min max])
  (:require [quil.core :as q]
            [quil.middleware :as middleware]
            [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :refer :all]))

(defn linspace [a b n]
  (let [step (/ (- b a) (dec n))]
    (range a (+ b step) step)))

(defn catmull-rom-spline
  ([p0 p1 p2 p3] (catmull-rom-spline p0 p1 p2 p3 20 0.5))
  ([p0 p1 p2 p3 n-points alpha]
   (letfn [(tj [ti pi pj]
             (let [[xi yi] pi
                   [xj yj] pj]
               (+ (m/pow (m/pow (+ (m/pow (- xj xi) 2) (m/pow (- yj yi) 2)) 0.5) alpha) ti)))
           (bc-prod [m1 m2]
             (if (== (last (m/shape m1)) 1)
               (bc-prod m2 m1)
               (let [m1-bc (m/broadcast m1 [(first (m/shape m2)), (last (m/shape m1))])
                     m2-bc (map #(m/broadcast (first %) [(last (m/shape m1))]) m2)]
                 (* m1-bc m2-bc))))]
     (let [t0 0
           t1 (tj t0 p0 p1)
           t2 (tj t1 p1 p2)
           t3 (tj t2 p2 p3)
           t (linspace t1 t2 n-points)
           t' (m/reshape t [(count t) 1])
           a1 (+ (bc-prod (/ (- t1 t') (- t1 t0)) p0) (bc-prod (/ (- t' t0) (- t1 t0)) p1))
           a2 (+ (bc-prod (/ (- t2 t') (- t2 t1)) p1) (bc-prod (/ (- t' t1) (- t2 t1)) p2))
           a3 (+ (bc-prod (/ (- t3 t') (- t3 t2)) p2) (bc-prod (/ (- t' t2) (- t3 t2)) p3))
           b1 (+ (bc-prod (/ (- t2 t') (- t2 t0)) a1) (bc-prod (/ (- t' t0) (- t2 t0)) a2))
           b2 (+ (bc-prod (/ (- t3 t') (- t3 t1)) a2) (bc-prod (/ (- t' t1) (- t3 t1)) a3))]
       (+ (bc-prod (/ (- t2 t') (- t2 t1)) b1) (bc-prod (/ (- t' t1) (- t2 t1)) b2))))))

(defn catmull-rom-chain [points]
  (mapcat #(catmull-rom-spline %1 %2 %3 %4) points (drop 1 points) (drop 2 points) (drop 3 points)))

(defn setup []
  (q/frame-rate 30)
  (catmull-rom-chain [[100 100] [100 200] [200 300] [300 300] [300 200] [300 100] [100 100] [100 200] [200 300]]))

(defn update-state [state]
  state)

(defn draw-state [state]
  (dorun (map #(q/line %1 %2) state (drop 1 state))))


(q/defsketch catmull-rom-quil
  :title "catmull-rom"
  :size [500 500]
  :setup setup
  :update update-state
  :draw draw-state
  :settings #(q/smooth 5)
  :features [:keep-on-top]
  :middleware [middleware/fun-mode])
