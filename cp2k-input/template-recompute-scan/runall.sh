for a in *just*; do cd $a; nf=$(ls f*xyz | wc -l | awk '{print $1}'); if [ $nf -ge 1 ]; then sbatch ../run.sh; echo $a; fi; cd ..;  done
