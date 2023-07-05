#!/bin/bash -l
#SBATCH --job-name=cp2k-md
#SBATCH --account=s1108
#SBATCH --time=06:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread

module load CP2K

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

ulimit -s unlimited

cp ../cp2k.inp .

mkdir finished

for f in f-*.xyz; do

tail -n +3 $f > init.xyz
head -n 2 $f | tail -n 1 | sed 's/\"/ /g' | awk '{print "ABC",$2,$6,$10}' > init.ABC

srun -n $SLURM_NTASKS cp2k.psmp -i cp2k.inp -o cp2k.out
wait
mv $f finished

lx=$(grep 'ABC' init.ABC | awk '{print $2}')
ly=$(grep 'ABC' init.ABC | awk '{print $3}')
lz=$(grep 'ABC' init.ABC | awk '{print $4}')

paste H2O-TiO2-pos-1.xyz H2O-TiO2-frc-1.xyz | awk -v lx=$lx -v ly=$ly -v lz=$lz 'NF==2{print $1} $1=="Ti" || $1=="O" || $1=="H" {print $1,$2,$3,$4,$6,$7,$8} $1=="i"{printf("Lattice=\"%4.4f 0 0 0 %4.4f 0 0 0 %4.4f\" energy=%4.8f Properties=species:S:1:pos:R:3:force:R:3\n", lx,ly,lz, $NF)}' > finished/cp2k-$f

done

wait
