(ql:quickload :gsll)
(ql:quickload :ieee-floats)
(asdf:load-system :cl-libplot)

(defparameter *repository-address* "http://xtof.disque.math.cnrs.fr/data/")
(defparameter *data-file-names* (let ((prefix "Locust_")
				      (suffix ".dat")
				      (numbers '("1" "2" "3" "4")))
				  (mapcar #'(lambda (i) (concatenate 'string prefix i suffix))
					  numbers)))
(export '*data-file-names*)
(defparameter *compressed-data-file-names* (mapcar #'(lambda (n) (concatenate 'string n ".gz"))
						   *data-file-names*))
(defparameter *full-data-file-names* (mapcar #'(lambda (n) (concatenate 'string *repository-address* n))
					     *compressed-data-file-names*))
(mapcar #'(lambda (x) (sb-ext:run-program "/usr/bin/wget" (list x))) *full-data-file-names*)
(mapcar #'probe-file *compressed-data-file-names*)
(mapcar #'(lambda (x) (sb-ext:run-program "/bin/gunzip" (list x))) *compressed-data-file-names*)
(mapcar #'probe-file *data-file-names*)

(defpackage "SPIKE-O-MATIC"
  (:use :cl :gsll :ieee-floats :cl-libplot)
  (:nicknames "SOM"))

(in-package "SOM")

(defclass mts ()
  ((data :accessor mts-data :initarg :data :documentation "The data matrix.")
   (sampling-rate :accessor mts-sr :initarg :sr :initform 15000 :documentation "The sampling rate (in s) used during the recording.")
   (start :accessor mts-start :initarg :start :initform 0 :documentation "The time (in s) at which the recording started.")
   (nb-channels :accessor mts-nbc :initarg :nbc :initform 4 :documentation "The number of channels to process simultaneoiusly.")
   (duration :accessor mts-duration :initarg :duration :documentation "The recording duration (in s)."))
  (:documentation "A mutivariate time-seris class containing the data in a matrix form (slot data) with one channel per column.
Contextual information like the sampling rate (slot sampling-rate), the time at which the recording started (slot start), 
the number of channels to process as a whole (slot nb-channels) and the recording duration (slot duration) is also included.
The time unit is the second."))

;;(import 'common-lisp-user::*data-file-names*)

(defun mk-mts-from-double (data-file-names-lst &key (sampling-rate 15000) (start 0))
  "Returns an mts object whose data slot contains the matrix made out of 
the content of the files named in the list data-file-names-lst"
  (let* ((nbc (length data-file-names-lst))
	 (data-length (with-open-file (in (nth 0 data-file-names-lst)
					  :direction :input
					  :element-type '(unsigned-byte 64))
			(file-length in)))
	 (data-dim (list data-length nbc))
	 (data (grid:make-foreign-array 'double-float :dimensions data-dim :initial-element 0d0)))
    (do ((ch-idx 0 (1+ ch-idx)))
	((> ch-idx (1- nbc)))
      (with-open-file (in (nth ch-idx data-file-names-lst)
			  :direction :input
			  :element-type '(unsigned-byte 64))
	(dotimes (i data-length)
          (grid:gsetf (grid:aref data i ch-idx) (ieee-floats:decode-float64 (read-byte in))))))
    (make-instance 'mts :data data :sr sampling-rate :start start :nbc nbc :duration (coerce (/ data-length sampling-rate) 'double-float))))

(defparameter *locust-data* (mk-mts-from-double cl-user:*data-file-names*))

(defun show (mts-obj &key (min nil) (max nil) (BITMAPSIZE "800x800+0+0") (line-width 0.01d0) (pen-color-name "black"))
  "Generates a 'bare-bone' plot of its first argument which must be a mts object.
Argument min should be the floor of the data slot of the first argument (nil is the default).
Argument max should be the ceiling of the data slot of the first argument (nil is the default).
Arguments BITMAPSIZE, line-width and pen-color-name are passed to the libplot plotter.
The function generates a plot as a side effect and returns a plist with elements:
  ':plotter-par', the pointer to the plotter parameter structure
  ':plotter', the pointer to the plotter structure
These two pointers have to be deleted by the user!"
  (if (not min) (setf min (floor (iter:iter (iter:for e :matrix-element (mts-data mts-obj)) (iter:minimize e)))))
  (if (not max) (setf max (ceiling (iter:iter (iter:for e :matrix-element (mts-data mts-obj)) (iter:maximize e)))))
  (let* ((range (coerce (- max min) 'double-float))
	 (plotter-params (cllp:newplparams))
	 (nbc (mts-nbc mts-obj))
	 (nbs (car (grid:dimensions (mts-data mts-obj)))))
    (cllp:setplparam plotter-params "BITMAPSIZE" BITMAPSIZE)
    (let ((plotter (cllp:newpl "X"
			       (cffi:foreign-symbol-pointer "stdin")
			       (cffi:foreign-symbol-pointer "stdout")
			       (cffi:foreign-symbol-pointer "stderr")
			       plotter-params))
	  (offset 0d0))
      (cllp:openpl plotter)
      (cllp:flinewidth plotter line-width)
      (cllp:pencolorname plotter pen-color-name)
      (cllp:lp-fspace plotter 0d0 0d0 (coerce nbs 'double-float) (* nbc range))
      (dotimes (ch-idx nbc)
	(setf offset (- (* (- (1- nbc) ch-idx) range) min))
	(cllp:fmove plotter 0d0 offset)
	(dotimes (i nbs)
	  (cllp:fcont plotter (coerce i 'double-float) (+ (grid:aref (mts-data mts-obj) i ch-idx) offset)))
	)
      (cllp:closepl plotter)
      (list :plotter-par plotter-params :plotter plotter))))

(defparameter *plot* (show *locust-data* :min -10 :max 12))

(cllp:deletepl (getf *plot* :plotter))
(cllp:deleteplparams (getf *plot* :plotter-par))

(defmethod quantiles ((data mts) (prob list))
  "Computes the quantiles of the data slot of its first argument
at each probability element of the list making up its second argument.
A grid array is returned with as many columns as columns in the original
mts object and as many rows as there are probabilities in the list."
  (mapcar #'(lambda (p) 
                (if (or (< p 0d0) (< 1.0d0 p))
                    (error "A probability p must satisfy 0 ≤ p ≤ 1!"))) 
            prob)
  (let* ((len (car (grid:dimensions (mts-data data))))
	 (nbc (mts-nbc data))
	 (res (make-array (list (length prob) nbc) :element-type 'double-float :initial-element 0d0))
	 (col (grid:make-foreign-array 'double-float :dimensions len :initial-element 0d0))
	 (quant prob))
    (dotimes (j nbc)
      (setf col (grid:column (mts-data data) j))
      (setf col (gsll:sort-vector col))
      (setf quant (if (= (length prob) 1)
		      (list (gsll:quantile col (car prob)))
		      (mapcar #'(lambda (p) (gsll:quantile col p)) prob)))
      (dotimes (i (length quant))
	(setf (aref res i j) (nth i quant))))
    res))

(quantiles *locust-data* '(0d0 0.25d0 0.5d0 0.75d0 1d0))

(iter:iter (iter:for col :matrix-column (mts-data *locust-data*)) (princ (standard-deviation col)) (princ " "))

(defun diff (x &key (lag 1))
    (declare (optimize (speed 3)))
    (declare (type grid:vector-double-float x))
    (declare (fixnum lag))
    (let* ((len (- (car (grid:dimensions x)) lag))
           (res (grid:make-foreign-array 'double-float :dimensions len :initial-element 0d0)))
      (declare (type grid:vector-double-float res))
      (declare (fixnum len))
      (dotimes (i len res)
        (grid:gsetf (grid:aref res i) (- (grid:aref x (+ i lag)) (grid:aref x i))))))

(defun unique (x)
    (declare (optimize (speed 3)))
    (declare (type grid:vector-double-float x))
    (let* ((len (car (grid:dimensions x)))
           (sorted-x (grid:slice x (list (list ':range 0 (1- len)))))
           (res (grid:make-foreign-array 'double-float :dimensions len :initial-element 0d0))
           (j 0)
           (v 0d0))
      (declare (type grid:vector-double-float sorted-x res))
      (declare (double-float v))
      (declare (fixnum len j))
      (setf sorted-x (gsll:sort-vector sorted-x))
      (setf v (grid:aref sorted-x 0))
      (grid:gsetf (grid:aref res 0) v)
      (do* ((i 1 (1+ i)))
           ((> i (1- len)))
        (setf v (grid:aref sorted-x i))
        (cond ((> v (grid:aref res j))
               (setf j (1+ j))
               (grid:gsetf (grid:aref res j) v))))
      (grid:slice res (list (list ':range 0 j)))))

(iter:iter (iter:for col :matrix-column (mts-data *locust-data*)) (iter:collect (reduce #'min (grid:copy-to (diff (unique col))))))
