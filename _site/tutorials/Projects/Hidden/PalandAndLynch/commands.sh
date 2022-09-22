python extractCDS.py sequence.gb Proteins CDS

mkdir CDS_aligned ; cd CDS ; for i in * ; do muscle -in $i -out ../CDS_aligned/${i} ; done ; cd ..

python /home/boussau/Work/scripts/CatAlignedFastaSeqsBEST.py CDS_aligned concatenate.fasta


