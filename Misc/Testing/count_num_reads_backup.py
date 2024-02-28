import pysam

in_dir = "/scratch/ProjectDir/aligns/"
files = ["DB181.cram",
         "DB182.cram",
         "DB183.cram",
         "DB184.cram",
         "DB185.cram",
         "DB186.cram",
         "DB187.cram",
         "DB188.cram",
         "DB189.cram",
         "DB190.cram",
         "DB191.cram",
         "DB192.cram",
         "DB193.cram",
         "DB194.cram",
         "DB195.cram",
         "DB196.cram",
         "DB197.cram",
         "DB198.cram",
         "DB199.cram",
         "DB200.cram",
         "DB201.cram",
         "DB202.cram",
         "DB203.cram",
         "DB204.cram",
         "DB205.cram",
         "DB206.cram",
         "DB207.cram",
]

total_files = len(files)
print(total_files)
ref_len = 3209286105

for cram_file in files:
    f = pysam.AlignmentFile(in_dir+cram_file)
    c = f.count()
    cov = (c*98)/ref_len
    print(str.split(cram_file, ".")[0], cov)

