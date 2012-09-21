(ql:quickload :gsll)
(ql:quickload :ieee-floats)
(asdf:clear-configuration)
(asdf:load-system :cl-libplot)

(load (compile-file "spike-sorting.lisp"))

(setf *print-length* 20)
(setf grid:*print-foreign-array-readably* nil)

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

(in-package "SOM")

(defparameter *locust-data* (mk-mts-from-double cl-user:*data-file-names*))

(defparameter *plot* (show *locust-data* :min -10 :max 12))

(quantiles *locust-data* '(0d0 0.25d0 0.5d0 0.75d0 1d0))

(iter:iter (iter:for col :matrix-column (mts-data *locust-data*)) (princ (standard-deviation col)) (princ " "))

(defparameter *step-sizes*
  (iter:iter (iter:for col :matrix-column (mts-data *locust-data*)) (iter:collect (reduce #'min (grid:copy-to (diff (unique col))))))) 

(defparameter *locust-data-i* (A-to-D-output *locust-data*))

(show *locust-data-i* :min -1352 :max 1584)

(show (window *locust-data-i* :from 0.2d0 :to 0.4d0) :min -1352 :max 1584)
