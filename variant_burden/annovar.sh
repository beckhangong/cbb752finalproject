# run Annovar to get Ensembl gene, dbsnp, and various other annotations
perl annovar/table_annovar.pl chr3.vcf annovar/humandb/ -vcfinput -buildver hg19 -out myanno -remove -nastring . --protocol ensGene,avsnp150,popfreq_max_20150413,clinvar_20190305 --operation g,f,f,f

# remove intergenic, intronic, and ncRNA-annotated variants, and also drop unneeded data columns from output (last several are VCF INFO fields)
awk '$6!="intergenic" && $6!~"ncRNA" && $6!="intronic"' ensembl_anno.hg19_multianno.txt | cut -f 1-18 >  chr3_annovar_out.txt

# after this, manually changed the name of last remaining data column from Otherinfo to variant_freq
