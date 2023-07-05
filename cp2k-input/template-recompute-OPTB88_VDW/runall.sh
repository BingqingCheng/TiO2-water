for a in *just*; do cd $a; if [ ! -e finished ]; then sbatch ../run.sh; echo $a; fi; cd ..;  done
