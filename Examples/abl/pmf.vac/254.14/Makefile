all: prepare cut 

cut: 
	mkdir -p history
	python cut.py 

sort: 
	python sort.py dist 0 15 
	python sort.py dist 1 15 
	python sort.py dist 2 15 
	python sort.py dist 3 15 
	python sort.py dist 4 15 
	python sort.py dist 5 15 
	python sort.py dist 6 15 
	python sort.py dist 7 15 
	python sort.py dist 8 15 
	python sort.py dist 9 15 
	python sort.py dist 10 15 
	python sort.py dist 11 15 
	python sort.py dist 12 15 
	python sort.py dist 13 15 
	python sort.py dist 14 15 
	python sort.py dist 15 15 
	python sort.py dist 16 15 
	python sort.py dist 17 15 
	python sort.py dist 18 15 
	python sort.py dist 19 15 

prepare: sort
	#sh prepare.com
	python prepare.py
