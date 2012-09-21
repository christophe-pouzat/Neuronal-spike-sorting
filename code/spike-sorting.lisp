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
	 (quant prob)
	 (is-integer (integerp (grid:aref (mts-data data) 0 0))))
    (dotimes (j nbc)
      (if is-integer
	  (dotimes (i len) (grid:gsetf (grid:aref col i) (coerce (grid:aref (mts-data data) i j) 'double-float)))
	  (setf col (grid:column (mts-data data) j)))
      (setf col (gsll:sort-vector col))
      (setf quant (if (= (length prob) 1)
		      (list (gsll:quantile col (car prob)))
		      (mapcar #'(lambda (p) (gsll:quantile col p)) prob)))
      (dotimes (i (length quant))
	(setf (aref res i j) (nth i quant))))
    res))

(defun show (mts-obj &key (min nil) (max nil) (BITMAPSIZE "800x800+0+0") (line-width 0.01d0) (pen-color-name "black"))
  "Generates a 'bare-bone' plot of its first argument which must be a mts object.
Argument min should be the floor of the data slot of the first argument (nil is the default).
Argument max should be the ceiling of the data slot of the first argument (nil is the default).
Arguments BITMAPSIZE, line-width and pen-color-name are passed to the libplot plotter.
The function generates a plot as a side effect and returns a list with two integers
which should both be 0 if the 'external' structures created to make the plot were
successufully deleted."
  (if (not min) (setf min (floor (iter:iter (iter:for e :matrix-element (mts-data mts-obj)) (iter:minimize e)))))
  (if (not max) (setf max (ceiling (iter:iter (iter:for e :matrix-element (mts-data mts-obj)) (iter:maximize e)))))
  (let* ((range (- max min))
	 (plotter-params (cllp:newplparams))
	 (nbc (mts-nbc mts-obj))
	 (nbs (car (grid:dimensions (mts-data mts-obj)))))
    (cllp:setplparam plotter-params "BITMAPSIZE" BITMAPSIZE)
    (cond ((floatp (grid:aref (mts-data mts-obj) 0 0))
	   (let ((plotter (cllp:newpl "X"
				      (cffi:foreign-symbol-pointer "stdin")
				      (cffi:foreign-symbol-pointer "stdout")
				      (cffi:foreign-symbol-pointer "stderr")
				      plotter-params))
		 (offset 0d0))
	     (setf range (coerce range 'double-float))
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
	     (list (cllp:deleteplparams plotter-params) (cllp:deletepl plotter)))) 
	((integerp (grid:aref (mts-data mts-obj) 0 0))
	 (let ((plotter (cllp:newpl "X"
				    (cffi:foreign-symbol-pointer "stdin")
				    (cffi:foreign-symbol-pointer "stdout")
				    (cffi:foreign-symbol-pointer "stderr")
				    plotter-params))
	       (offset 0))
	   (cllp:openpl plotter)
	   (cllp:flinewidth plotter line-width)
	   (cllp:pencolorname plotter pen-color-name)
	   (cllp:lp-space plotter 0 0 nbs (* nbc range))
	   (dotimes (ch-idx nbc)
	     (setf offset (- (* (- (1- nbc) ch-idx) range) min))
	     (cllp:move plotter 0 offset)
	     (dotimes (i nbs)
	       (cllp:cont plotter i (+ (grid:aref (mts-data mts-obj) i ch-idx) offset)))
	     )
	   (cllp:closepl plotter)
	   (list (cllp:deleteplparams plotter-params) (cllp:deletepl plotter))))
	(t (cllp:deleteplparams plotter-params)))))

(defun diff (x &key (lag 1))
  "Returns x[i]-x[i-lag] where 'x' is the first argument which must be
a 'grid double float vector' and 'lag' is a positive integer. A 
'grid double float vector' is returned."
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
  "Returns a 'grid:vector-double-float' containing only once each value present
in its argument which must also be a 'grid:vector-double-float'."
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

(defun A-to-D-output (mts-obj)
  "Returns a MTS object identical to the one given as argument, except that
each column of the 'data' slot has been divided by the discretization step
amplitude and rounded to become a 16 bit integer."
  (if (not (eq 'MTS (type-of mts-obj))) ;; check argument
      (error "The argument should be a MTS object."))
  (let ((step-sizes (iter:iter (iter:for col :matrix-column (mts-data mts-obj))
			       (iter:collect (reduce #'min (grid:copy-to (diff (unique col))))))) ;; discretization steps
	(nbc (mts-nbc mts-obj)) ;; numbre of channels
	(len (car (grid:dimensions (mts-data mts-obj)))) ;; number of sampling points per channel
	(mat (grid:make-foreign-array '(SIGNED-BYTE 16) :dimensions (grid:dimensions (mts-data mts-obj)) :initial-element 0)))
    (dotimes (j nbc) ;; loop on columns
      (let ((divisor (nth j step-sizes))) ;; get discretization step of the column
	(dotimes (i len) (grid:gsetf (grid:aref mat i j) (round (/ (grid:aref (mts-data mts-obj) i j) divisor))))))
    (make-instance 'mts :data mat :sr (mts-sr mts-obj) :duration (mts-duration mts-obj) :start (mts-start mts-obj) :nbc nbc)))

(defun window (mts-obj &key (from 0d0) (to 0.2d0))
  "Returns a MTS object obtained by selecting a time window between arguments
'from' and 'to' (in s) of its first argument. Inspired by R window function."
  (if (not (eq 'MTS (type-of mts-obj))) ;; check argument
      (error "The argument should be a MTS object."))
    (let* ((sr (mts-sr mts-obj))
	   (len (car (grid:dimensions (mts-data mts-obj))))
	   (from-d (max 0 (round (* from sr))))
	   (to-d (min (1- len) (round (* to sr))))
	   (idx-select (list (list ':range from-d (1- to-d)) ':all))
	   (mat (grid:slice (mts-data mts-obj) idx-select)))
      (make-instance 'mts :data mat :sr sr :duration (/ (- to-d from-d) sr) :start (/ from-d sr) :nbc (mts-nbc mts-obj))))
