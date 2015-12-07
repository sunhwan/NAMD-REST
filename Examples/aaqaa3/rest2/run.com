qsub -t 1:00:00 -n 128 -A US-REST2 --env icnt=0:n=8:nrep=16 --mode script run.sh
