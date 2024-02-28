#!/bin/bash
## java install
sudo apt-get update
sudo apt-get upgrade -y
sudo apt-get install openjdk-11-jdk
sudo apt-get install -y libfuse-dev fuse curl wget software-properties-common

## AWS Mounting Drive
sudo mkfs.ext4 /dev/nvme0n1
sudo mkdir /mnt/nvme0n1
sudo mount /dev/nvme0n1 /mnt/nvme0n1
sudo chmod 777 /mnt/nvme0n1

mkdir -p /home/ubuntu/.local/bin

echo "" >> .bashrc
echo "export PATH=$PATH:$HOME/.local/bin" >> $HOME/.bashrc
echo "export PYTHONPATH=$PYTHONPATH:$HOME/.local/bin" >> $HOME/.bashrc
echo "export STORAGE_PROFILE=collab" >> $HOME/.bashrc
source .bashrc

## score client
wget -O score-client.tar.gz https://artifacts.oicr.on.ca/artifactory/dcc-release/bio/overture/score-client/\[RELEASE\]/score-client-\[RELEASE\]-dist.tar.gz
tar -xvzf score-client.tar.gz
# add key to score client
rm score-client*/conf/application.properties
mv application.properties score-client*/conf/
scli="/home/ubuntu/score-client*/bin/score-client"
export STORAGE_PROFILE=collab

## install samtools
sudo apt install samtools

## install vim
sudo apt install vim

## install sshpass for scp with password send
sudo apt install sshpass

### optional install docker
## sudo apt install docker
## sudo docker pull overture/score
## sudo docker run -it overture/score
## sudo apt-get update
## sudo apt-get upgrade


## install dependancies for tools
sudo apt-get install python3-pip
pip3 install Cython numpy pysam
pip3 install click
pip3 install click>=7.1.2 -U
pip3 install numpy>1.18.5 -U

sudo apt install autoconf
sudo apt install libbz2-dev liblzma-dev

sudo apt install minimap2

#export PYTHONPATH=$PYTHONPATH:$HOME/.local/bin
#export PATH=$PATH:$HOME/.local/bin

## install dysgu # clone for coverage script
pip install dysgu
## build from source 
# cd /home/ubuntu/
git clone --recursive https://github.com/kcleal/dysgu.git
# cd dysgu/dysgu/htslib
# autoreconf -i
# ./configure
# make
# sudo make install
# cd ../../
# pip install -r requirements.txt
# pip install .
# bash INSTALL.sh # errors
#export PATH=$PATH:$HOME/.local/bin

mkdir -p /mnt/nvme0n1/dysgu_out
sudo chmod 777 /mnt/nvme0n1/dysgu_out
cd /home/ubuntu/

## install teltool
git clone https://github.com/kearseya/teltool.git
cd teltool
pip3 install -e .
cd /home/ubuntu/
mkdir -p /mnt/nvme0n1/teltool_out


## download refernce
mkdir -p /mnt/nvme0n1/reference_genomes
cd /mnt/nvme0n1/reference_genomes
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz

## download data (terraform touch ${s}.sample method)
#cd /home/ubuntu/
#sample_id="*.sample"
#${scli} download --object-id $(basename ${sample_id} .sample) --output-dir /home/ubuntu


## run dysgu 
mkdir /mnt/nvme0n1/dysgu_out
ref="/mnt/nvme0n1/reference_genomes/hs37d5.fa"
cov="python3 /home/ubuntu/dysgu/scripts/coverage2bed.py"
readarray -t s < /home/ubuntu/nomini_object_ids


## from previous mounting method
#s=(`ls *.bam`)
#s=$(cat filter/bampaths | tr "\n" " ")

echo -e "--------------------\n\tLIST (${#s[@]})\n\n${s[*]}\n"
echo -e "--------------------\n\n\tSTART\n\n--------------------"
for s in ${s[*]}
	do
	## download
	cd /mnt/nvme0n1/
	${scli} download --object-id ${s} --output-dir /mnt/nvme0n1
	## run dysgu
	i="/mnt/nvme0n1/*.bam"
	b=$(basename ${i} .bam)
	echo "Running dysgu ${b}"
	## run
	/home/ubuntu/.local/bin/dysgu run -v2 --metrics --low-mem --exclude /home/ubuntu/hg19-blacklist.v2.bed --mq 15 --max-cov -1 -x --clip-length 30 --min-support 5 --search /home/ubuntu/chromosome.bed ${ref} tmp_${b} ${i} -o ${b}.vcf 2> ${b}.log ; ${cov} --out-bin-size 1000 -w tmp_${b} > ${b}_cov.bed ; rm -rf tmp_${b}
	## fetch call
	#/home/ubuntu/.local/bin/dysgu fetch -p2 --exclude /home/ubuntu/hg19-blacklist.v2.bed --mq 15 --max-cov -1 -x --clip-length 30 --search /home/ubuntu/chromosome.bed tmp_${b} ${i} 2> ${b}.log ; /home/ubuntu/.local/bin/dysgu call -p 1 -v2 --metrics --mq 15 --max-cov -1 -x --clip-length 30 --min-support 5 --low-mem --ibam ${i} ${ref} tmp_${b}/*.dysgu_reads.bam -o ${b}.vcf 2>> ${b}.log; ${cov} --out-bin-size 1000 -w tmp_${b} > ${b}_cov.bed ; rm -rf tmp_${b}

	## move and transfer output
	mv ${b}.vcf /mnt/nvme0n1/dysgu_out/
	mv ${b}.log /mnt/nvme0n1/dysgu_out/
	mv ${b}_cov.bed /mnt/nvme0n1/dysgu_out/
	sshpass -f "/home/ubuntu/pw.txt" scp -r /mnt/nvme0n1/dysgu_out/* User@Host:/scratch/User/icgc/dysgu
	## run teltool
	if grep -Fxq ${b} noteltool
	then
		echo "Trimming ${b}"
		teltool trim -i ${i} -o teltool_out
		mv coverages.csv teltool_out/${b}_cov.csv
		mv *_tel.bam teltool_out/${b}_tel.bam
		mv *_tel.bam.bai teltool_out/${b}_tel.bam.bai
		sshpass -f "/home/ubuntu/pw.txt" scp -r /mnt/nvme0n1/teltool_out/* User@Host:/scratch/User/icgc/teltool
	else
		echo "Already done teltool ${b}"
	fi
	
	## remove sample for space
	rm *.bam*
## run jobs in parallel
#	echo ${b}
#	jobs=($(jobs -p))
#	while (( ${#jobs[*]} >= 2 ))
#	do
#		sleep 30
#		jobs=($(jobs -p))
#	done
done

