


# -------
# There is a split csv and the gc is sent combined with full gc df
# it then queries the df, makes a fasta, runs gctree, and finally does a second wrangle
# THIS IS IF GC != 0
# The other case will 

GCDF=gc-df.csv                  # Full PR cell HK database (csv)
GC_NUM_HC=87                    # GC you would like to extract
QUERY="GC_num_HC==${GC_NUM_HC}" # The query on the dataset before fasta extraction

# query the gc dataframe
python gcreplay-tools.py query-df -df $GCDF -q $QUERY -o gc_${GC_NUM_HC}_hk.csv
echo done with query

# convert some columns into headers and sequences for a basic fasta format conversion
python gcreplay-tools.py gc-df-to-fasta --gc-df gc_${GC_NUM_HC}_hk.csv -o gc_${GC_NUM_HC}_hk.fasta


# -------
# Now running gctree
# note this is not including the HDAG stuff yet

#MUTFILES_1=../../data/mutability-files/
#MUTFILES_2=../../data/mutability-files/
# TODO DEDUP

#-> CWD
#     |-> gc13HK
#     |     |-> gc13HK.fa
#     |-> gc18HK
#     |     |-> gc18HK.fa
#     |-> gc21HK
#     |     |-> gc21HK.fa
#   ...
#
## 
#cd gc[N]HK

## TODO describe
## First we deduplicate the sequences and convert to phylip alignment format, 
## and also determine sequence abundances. The deduplicate command writes 
## the phylip alignment of unique sequences to stdout (which we redirect to a file here). 
#deduplicate gc_${GC_NUM_HC}_hk.fasta \
#    --root naive \
#    --idmapfile idmap-gc${GC_NUM_HC}-HK.txt \ # HERE What is the idmap. ID_HK -> ?
#    --abundance_file abundances.csv > deduplicated.phylip
#
## TODO describe
#mkconfig deduplicated.phylip dnapars > dnapars.cfg
#
## TODO describe
#dnapars < dnapars.cfg > dnapars.log
#
## TODO describe
#gctree infer outfile abundances.csv \
#    --root naive \
#    --frame 1 \
#    --frame2 1 \
#    --positionmapfile ../210929PR_H_index-imgtTs.txt \ # position map file
#    --positionmapfile2 ../210930PR_K_index-imgtTs.txt \
#    --chain_split 336 \
#    --idlabel \
#    --verbose \
#    | tee gctree.inference.log

# TODO describe
## Post-GCtree processing
#- import idmap, add to GCdf
#- import putative nodes, fill idmap info, add to GCdf
#- calculate total mut#
#- calculate DMS scores
