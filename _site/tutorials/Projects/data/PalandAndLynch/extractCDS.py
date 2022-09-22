from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import collections
import os
import sys


gb_file = sys.argv[1]
output_protein_folder = sys.argv[2]
output_cds_folder = sys.argv[3]

os.mkdir(output_protein_folder)
os.mkdir(output_cds_folder)

gene = ""
sp = ""
prot_id = ""
aa_seq = ""
first_time = True
ignore = False

with open(gb_file,  "r") as handle:
    for record in SeqIO.parse(handle, "genbank"):
        print(record)
        sp = str(record.annotations['organism'])
        strain = ""
        description = str(record.description)
        strain = description.split()[3]
        print(strain)
        sp_strain = sp.replace(" ", "_") + "_" + strain
        all_aas = collections.defaultdict(int)
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    try:
                        print(feature)
                        #prot_id = str(feature.qualifiers["product"])
                        print(feature.qualifiers["gene"][0])
                        gene = str(feature.qualifiers["gene"][0])
    #                    print(feature.qualifiers["protein_id"])
                        aa_seq = str(feature.qualifiers["translation"][0])
                        cds_seq = str(feature.extract(record).seq)
                        print(cds_seq)
                        #print(">" + sp + " " +gene + " " + prot_id + "\n" + aa_seq)
                        with open(os.path.join(output_protein_folder, gene+".fasta"), "a") as output_handle:
                            output_handle.write(">" + sp_strain + "\n" + aa_seq + "\n")
                        with open(os.path.join(output_cds_folder, gene+".fasta"), "a") as output_handle:
                            output_handle.write(">" + sp_strain + "\n" + cds_seq + "\n")
                        #print(all_aas)
                    except:
                        print("PROBLEM with sp "+ sp )
                        #print ( "\n\n and record: " + record)
                        print (" IGNORING")
                        ignore = True
                else:
                    pass
                    #print(feature)
