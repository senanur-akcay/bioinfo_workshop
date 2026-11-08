#analiz ortamını oluşturma
conda deactivate #rnaseq ortmından çıktım
conda env create -f envs/variants.yml
conda env list
#toolların versiyonlarını kontrol et
conda activate variants
bwa 2>&1 | grep Version
samtools --version | head -1
bcftools --version | head -1
which bgzip tabix

#data indirme
bash setup/day2.sh
#data kontrol
ls -lh data/

#okuma satır kontrolü
zcat data/NA12878_window_1.fastq.gz | head -4
for r in 1 2; do
  echo "R$r: $(( $(zcat data/NA12878_window_$r.fastq.gz | wc -l) / 4 )) reads"
done
zcat data/NA12878_window_1.fastq.gz | head -40000 | awk 'NR % 4 == 2 { print length($0) }' | sort -n | uniq -c | sort -k1,1nr | head -3

#FASTQC
mkdir -p results/fastqc
fastqc -t 4 --extract -o results/fastqc data/NA12878_window_1.fastq.gz data/NA12878_window_2.fastq.gz

#örnek sonuç yazdırma, FASTQC
cat results/fastqc/*_fastqc/summary.txt | sort -k2,2 -t$'\t' | cut -f1,2,3


#FASTp, dosyada değişiklik yapan ilk QC adımı #bu adımda sadece adaptörler kırpıldı, asıl kırpma sonrasında gerçekleştirilecektir.
mkdir -p results/fastp
fastp -i data/NA12878_window_1.fastq.gz -I data/NA12878_window_2.fastq.gz \
      -o results/fastp/NA12878_1.trim.fastq.gz -O results/fastp/NA12878_2.trim.fastq.gz \
      --detect_adapter_for_pe -w 4 \
      -j results/fastp/NA12878.fastp.json -h results/fastp/NA12878.fastp.html

#tekrar FASTQC
mkdir -p results/fastqc_trimmed
fastqc -t 4 --extract -o results/fastqc_trimmed results/fastp/NA12878_1.trim.fastq.gz results/fastp/NA12878_2.trim.fastq.gz

#HANGİLERİ GEÇTİ/HANGİLERİ KALDI
cat results/fastqc_trimmed/*_fastqc/summary.txt | sort -k2,2 -t$'\t' | cut -f1,2,3

#multiqc
multiqc -q -o results/multiqc results/fastqc results/fastqc_trimmed results/fastp

#referans ve index indirme (index referansta aradığın kromozom parçasının nerede old. gösterir.)
mkdir -p results/ref
gunzip -c data/chr20.fa.gz > results/ref/chr20.fa
samtools faidx results/ref/chr20.fa
cut -f1,2 results/ref/chr20.fa.fai
bwa index results/ref/chr20.fa 2>&1 | tail -3
ls -lh results/ref/

#BWA ile hizalama
mkdir -p results/align
bwa mem -t 4 -R '@RG\tID:NA12878\tSM:NA12878\tPL:ILLUMINA' results/ref/chr20.fa \
    results/fastp/NA12878_1.trim.fastq.gz results/fastp/NA12878_2.trim.fastq.gz 2> results/align/bwa.log \
  | samtools fixmate -m - - \
  | samtools sort -@ 4 -o results/align/NA12878.sorted.bam -   #binary alignment map (bam) elimdeki veri referansta var mı, nerede?
tail -3 results/align/bwa.log
ls -lh results/align/

#duplikatlar/indexler inceleme
samtools markdup -s results/align/NA12878.sorted.bam results/align/NA12878.bam 2>&1 | head -20
samtools index results/align/NA12878.bam

samtools flagstat results/align/NA12878.bam

#BAN dosyasının birkaç satırı
samtools view -H results/align/NA12878.bam | grep -v '^@SQ' | head
samtools view -H results/align/NA12878.bam | grep -c '^@SQ'
samtools view results/align/NA12878.bam chr20:10500000-10500300 | head -3 | cut -f1-9

#kapsam ve derinlik
samtools coverage -r chr20:10000000-12000000 results/align/NA12878.bam
samtools view -F 0xB04 results/align/NA12878.bam | cut -f5 > results/align/mapq.txt
sort -n results/align/mapq.txt | uniq -c | sort -k1,1nr | head -5

#Önemli not: kullandığınız panel ve sekanslama ne olursa olsun whole genome'a allign yapma daha uygun olur. false allignment önlenir.

#elde edilen hizlanmış ve sıralanmış dosyaları varyant bulma için hizalama
mkdir -p results/variants #mpileup varyant olasılığını hesaplar. 
bcftools mpileup -f results/ref/chr20.fa -r chr20:10000000-12000000 -a FORMAT/AD,FORMAT/DP -Ou results/align/NA12878.bam \
  | bcftools call -m -v -a GQ -Oz -o results/variants/NA12878.raw.vcf.gz
bcftools index -t results/variants/NA12878.raw.vcf.gz

#VCF dosyasını okumak
bcftools view -h results/variants/NA12878.raw.vcf.gz | grep -E '^##(INFO|FORMAT)=<ID=(DP|MQ|GT|GQ|AD|PL),'
bcftools view -h results/variants/NA12878.raw.vcf.gz | tail -1
bcftools head -h 0 -n 3 results/variants/NA12878.raw.vcf.gz | cut -f1-10

#Ne kadar SNP/Indel var? #henüz filtrelenmedi, filtrelemeden önce grafikte bakılır.
echo "all:    $(bcftools view -H results/variants/NA12878.raw.vcf.gz | wc -l)"
echo "SNPs:   $(bcftools view -H -v snps results/variants/NA12878.raw.vcf.gz | wc -l)"
echo "indels: $(bcftools view -H -v indels results/variants/NA12878.raw.vcf.gz | wc -l)"

#raw_variants dosyasini olusturacak
bcftools norm -f results/ref/chr20.fa -m -any -Oz -o results/variants/NA12878.raw.norm.vcf.gz results/variants/NA12878.raw.vcf.gz
bcftools index -t results/variants/NA12878.raw.norm.vcf.gz
bcftools query -H -f '%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%INFO/DP\t%INFO/MQ[\t%GT\t%GQ\t%AD]\n' \
  results/variants/NA12878.raw.norm.vcf.gz > results/variants/raw_calls.tsv
head -4 results/variants/raw_calls.tsv

#Filtreleme asamalari
bcftools filter -s LowQual -e 'QUAL<30 || INFO/DP<10 || INFO/DP>100 || INFO/MQ<30' \
  -Oz -o results/variants/NA12878.filtered.vcf.gz results/variants/NA12878.raw.norm.vcf.gz
bcftools index -t results/variants/NA12878.filtered.vcf.gz
bcftools view -H results/variants/NA12878.filtered.vcf.gz | cut -f7 | sort | uniq -c
bcftools view -f PASS -Oz -o results/variants/NA12878.pass.vcf.gz results/variants/NA12878.filtered.vcf.gz
bcftools index -t results/variants/NA12878.pass.vcf.gz

#filtrelemeyi gecenlerin istatsitik degerlerini dosyaya yazdir
bcftools stats results/variants/NA12878.pass.vcf.gz > results/variants/pass.stats.txt
grep -E '^SN' results/variants/pass.stats.txt | cut -f3-
grep -E '^TSTV' results/variants/pass.stats.txt | cut -f3-5

#Filtrelenmis olanlari tek bir dosyaya topla
bcftools query -f '%POS\t%QUAL\t%INFO/DP\t%FILTER\n' results/variants/NA12878.filtered.vcf.gz > results/variants/filtered_calls.tsv
head -3 results/variants/filtered_calls.tsv
#bu noktada r scripte geçip grafik oluştur.





















