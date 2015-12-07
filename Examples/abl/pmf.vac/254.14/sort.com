#~/local/namd/2.9/sortreplicas output/%s/dist.job0 48 1 2000000
#~/local/namd/2.9/sortreplicas output/%s/dist.job1 48 1 3000000
python sort.py dist 0 32
python sort.py dist 1 32
python sort.py dist 2 32
python sort.py dist 3 32
python sort.py dist 4 32
python sort.py dist 5 32
python sort.py dist 6 32
python sort.py dist 7 32
python sort.py dist 8 32
python sort.py dist 9 32
python prepare.py
