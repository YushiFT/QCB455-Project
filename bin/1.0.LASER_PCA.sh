#!/bin/bash
#SBATCH --job-name=laserPCA                # create a short name for your job
#SBATCH --nodes=1                          # node count
#SBATCH --ntasks-per-node=1                # total number of tasks across all nodes
#SBATCH --cpus-per-task=16                 # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem=60G                          # memory per node
#SBATCH --time=48:00:00 --qos=1wk          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin                  # send mail when process begins
#SBATCH --mail-type=end                    # send email when job ends
#SBATCH --mail-user=yushit@princeton.edu

# Extract PCs from genotype data of 1000 Genome Project
# define input
in_vcf=/scratch/tmp/yushi/QCB455/data/GEUVADIS_gen.vcf
in_laser=/Genomics/grid/users/yushit/.local/bin/LASER-2.04
in_vcf2geno=/Genomics/grid/users/yushit/.local/bin/LASER-2.04/vcf2geno
# define output
out_geno=/scratch/tmp/yushi/QCB455/temp/GEUVADIS_gen.conf
out_geno_pc_csv=/scratch/tmp/yushi/QCB455/data/GEUVADIS_gen_pcs.csv

# convert data format
echo 'converting VCF into CONF for LASER...'
$in_vcf2geno/vcf2geno --inVcf $in_vcf --out $out_geno

# extracting top 375 principal components
echo 'extracting top PCs...'
$in_laser/laser -p $out_geno -pca 1 -k 375 -o $out_geno_pc_csv

rm -f in_vcf
rm -f out_geno
