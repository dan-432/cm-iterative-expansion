#!/bin/bash
#SBATCH --partition=aoraki
#SBATCH --job-name=zip_results_and_data_ie
#SBATCH --output=zip_results_and_data_%j.log
#SBATCH --error=zip_results_and_data_%j.err
#SBATCH --time=12:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB                  # Request more memory
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dacda462@student.otago.ac.nz

tar -czvf expansion_results.tar.gz results/ ; tar -czvf expansion_data.tar.gz data/