import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import re
from glob import glob
from tqdm import tqdm
from os.path import getsize


### resfinder
for f in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/resfinder/*/ResFinder_results_tab.txt'):
    _df = pd.read_csv(f,sep='\t',comment='#')
    if _df.shape[0]!=0:
        print(_df)
        break

#! No resistence genes found
###    RGI and CARD

genome2rgi = {}
for otab in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/rgi/*.txt'):
    name = otab.split('/')[-1].rsplit('.')[0]
    _df = pd.read_csv(otab,sep='\t')
    genome2rgi[name] = _df
info_df = pd.read_csv('/mnt/maple/thliao/data/protein_db/ARG_related/CARD_3.2.6/aro_categories_index.tsv',sep='\t',index_col=0)
HEADER = ['qaccver', 'saccver', 'pident', 'length', 'mismatch',
          'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
odir = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/blast/CARD'
genome2subdf = {}
for blastout in tqdm(glob(f'{odir}/*.blastout')):
    if getsize(blastout) == 0:
        continue
    df = pd.read_csv(blastout, sep='\t', header=None)
    df.columns = HEADER
    df = df.loc[df['evalue'] < 1e-3, :]
    subdf = df.sort_values('evalue', ascending=True).groupby(
        'qaccver').head(n=1)
    subdf.loc[:,'pid'] = [_.split('|')[1] for _ in subdf['saccver']]
    subdf = subdf.loc[subdf['pid'].isin(info_df.index),:]
    genome2subdf[blastout.split('/')[-1].rsplit('.')[0]] = subdf

### antismash
indir = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/antismash'
from glob import glob
from Bio import SeqIO
import re
from collections import defaultdict
def get_region(idir):
    region2info = defaultdict(dict)
    for rfile in glob(join(idir,'*.region*.gbk')):
        for record in SeqIO.parse(rfile,'genbank'):
            f = [f for f in record.features if f.type =='protocluster'][0]
            category = ';'.join(f.qualifiers['category'])
            location = f.qualifiers['core_location'][0]
            location = re.findall(r'\[(\d+):(\d+)\]',location)
            product = ';'.join(f.qualifiers['product'])

            all_l = ';'.join([_.qualifiers['locus_tag'][0] for _ in record.features if _.type=='CDS'])
            region2info[rfile.split('.')[-2]]['contig'] = record.id
            region2info[rfile.split('.')[-2]]['location'] = location
            region2info[rfile.split('.')[-2]]['product'] = product
            region2info[rfile.split('.')[-2]]['category'] = category
            region2info[rfile.split('.')[-2]]['locus_tag'] = all_l
    return region2info

genome2results = {}
for fdir in glob(f'{indir}/*'):
    genome2results[basename(fdir).split('_')[0]] = get_region(fdir)
locus2region_product = {}
for genome,r2d in genome2results.items():
    for r,d in r2d.items():
        l_all = d['locus_tag'].split(';')
        product = d['product']
        for l in l_all:
            locus2region_product[l]=product
ofile = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/antismash_merged.tsv'
with open(ofile,'w') as f1:
    for l,pro in locus2region_product.items():
        f1.write(f"{l}\t{pro}\n")
products = set([v['product'] for _ in genome2results.values() for v in _.values()])
d = pd.DataFrame(index=products,columns=set(genome2results))
for gid,v in genome2results.items():
    c = defaultdict(list)
    for region,_d in v.items():
        names = []
        for s,e in _d['location']:
            names.append(_d['contig'] + f':{s}-{e}') 
        c[_d['product']].append(';'.join(names))
    for product, n in c.items():
        d.loc[product,gid] = ';'.join(n)
ofile = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/antismash_merged_g2product.tsv'
with open(ofile,'w') as f1:
    for l,pro in locus2region_product.items():
        f1.write(f"{l}\t{pro}\n")
genome2results
### BacARscan 
gid2arscan = {}
for otab in glob("/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/BacARscan/*.tab"):
    gid = otab.split('/')[-1].split('.')[0]
    _df = parse_hmm_otab(otab)
    _df = _df.loc[_df[2]<=1e-20,:]
    if _df is None:
        continue
    gid2arscan[gid] = _df
    
### KEGG

from Bio import SeqIO
ko2g_df = pd.read_csv(f"/mnt/maple/thliao/data/protein_db/kegg/ko_info.tab",sep='\t',header=None,index_col=0)
ko2g = {ko.split(':')[-1]:str(v).split(';')[0].strip() for ko,v in ko2g_df[1].to_dict().items()}
ko2info = {ko.split(':')[-1]:i for ko,i in ko2g_df[1].to_dict().items()}
basedir = '/home-user/thliao/project/coral_ruegeria/nanopore_processing/'
from api_tools.tk import *
from glob import glob
import pandas as pd
def parse_o(inf):
    l2ko = {}
    for row in open(inf).read().strip().split('\n'):
        rows = row.split('\t')
        l2ko[rows[0]] = ';'.join(sorted(rows[1:]))
    return l2ko
dfs = []
for kofamout in tqdm(glob(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KOFAMSCAN/*.kofamout")):
    gid2kegg2locus_info = defaultdict(lambda :defaultdict(list))
    l2ko = parse_o(kofamout)
    for locus,ko_l in l2ko.items():
        genome = kofamout.split('/')[-1].rsplit('.')[0]
        for ko in ko_l.split(';'):
            gid2kegg2locus_info[genome][ko].append(locus)
    gid2kegg2locus_info = {genome:{ko:','.join(list(set(l_list))) for ko,l_list in _d.items()} for genome,_d in gid2kegg2locus_info.items()}
    sub_df = pd.DataFrame.from_dict(gid2kegg2locus_info, orient='index')
    sub_df.to_csv(kofamout.replace('.kofamout','_anno.tab'),sep='\t',index=1)
    dfs.append(sub_df)
final_df = pd.concat(dfs,axis=0).fillna('NA')
final_df.to_csv(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno.tab", sep='\t', index=1)

comp2name = {'GNM004006175':"Parasedimentitalea marina", 
             'GNM003443535':"R. sp. AD91A", 
             'GNM000011965':"R. pomeroyi DSS-3", 
             'GNM000014065':"R. sp. TM1040"}
query_faa = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/merged.faa'
os.system(f"cat {basedir}/canu_o/*/09_prokka/*.faa > {query_faa}")
for _ in comp2name:
    os.system(f"cat /mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/prokka_o/{_}/{_}.faa >> {query_faa}")
os.system(f"cd /mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ && makeblastdb -in merged.faa  -dbtype prot")
length_g = len(glob(f"{basedir}/canu_o/*/09_prokka/*.faa"))+4
from api_tools.tk import read_hmmsearch_tbl
cmds = []
K = ['K00031','K00030']+'''K17222+K17223+K17224+K17225+K22622+K17226+K17227'''.split('+')
prefix=f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/extra_KO/{length_g}G_'
for ko in K:
    cmd = f"hmmsearch -o /dev/null --tblout {prefix}{ko} -T 0 --cpu 5 /home-user/thliao/db/protein_db/kegg/v20230301/profiles/{ko}.hmm {query_faa}"
    if not exists(f"{prefix}{ko}"):
        cmds.append(cmd)
        

from bin.multiple_sbatch import sbatch_all        
sbatch_all(cmds,thread_per_tasks=5,prefix_name='hmm')
kegg_df = pd.read_csv(
    '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno.tab', sep='\t', index_col=0,low_memory=False)
locus2ko = {locus: (ko2g.get(ko, ko), ko2info.get(ko, ko), ko)
            for ko, _d in kegg_df.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}
for extra in K:
    ofile = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/extra_KO/{length_g}G_{extra}"
    _df = read_hmmsearch_tbl(ofile,apply_preset_filter=False)
    sub_df = _df.loc[(~_df['target name'].isin(locus2ko)) & (_df.loc[:,"E-value(full)"] <= 1e-10),:]
    gid2kos = defaultdict(list)
    for l in sub_df['target name'].unique():
        gid2kos[l.split('_')[0]].append(l)
    gid2ko = {gid:','.join(kos) for gid,kos in gid2kos.items()}
    if extra not in kegg_df.columns:
        match_l = [gid2ko.get(gid,'') for gid in kegg_df.index]
    else:
        existed_matchs = list(kegg_df[extra])
        match_l = [','.join(list(set(gid2kos.get(gid,[]) + [existed_m]))) if str(existed_m)!='nan' else gid2ko.get(gid,'')
                   for existed_m, gid in zip(existed_matchs,kegg_df.index)]
    kegg_df.loc[:,extra] = match_l
kegg_df.to_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index=1)



## IS (insertion sequences)  see single_command.py
# cmds = []
# fna_list = glob(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/rotated_seqs/*.fna")
# for fna in fna_list:
#     sid = realpath(fna.split('/'))[-3]
#     odir = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/rotated_seqs/{sid}"
#     cmd = f"python /home-user/thliao/software/ISEScan-1.7.2.3/isescan.py --seqfile {fna} --output {odir} --removeShortIS --nthread 20"
#     if not exists(f"{odir}/{sid}.fna.csv"):
#         cmds.append(cmd)

genome2gbk = {}
for gbk in tqdm(glob(f"{basedir}/canu_o/*/09_prokka/*.gbk")):
    sid = gbk.split('/')[-3]
    if sid == gbk.split('/')[-1].replace('.gbk',''):
        genome2gbk[sid] = {_.id:_ for _ in SeqIO.parse(gbk, 'genbank')}
        
anno_df = pd.read_csv(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv",sep='\t',index_col=0,low_memory=True)
gids = [_.split('/')[-1].split('.')[0] for _ in genome2gbk]
sub_annodf = anno_df.reindex(gids)
locus2ko = {locus: (ko2g.get(ko,ko),
                    ko2info.get(ko,ko),
                    ko)
            for ko, _d in sub_annodf.to_dict().items()
            for genome, locus_list in _d.items()
        for locus in str(locus_list).split(',')}

MC_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = MC_df['MC'].to_dict()

IS_info_df = pd.DataFrame(columns=['cluster','start','stop','strand',
                                   'family','length of spanned region','contained genes','annotated kos','kos'])
IS_ID2locus_list = {}
for gff in tqdm(glob(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ISEScan/*/*/*.gff")):
    gff_df = pd.read_csv(gff,comment='#',sep='\t',header=None)
    s = gff.split('/')[-1].replace('.fna.gff','')
    if s not in genome2gbk:continue
    contig2seq = genome2gbk[s]
    sub_df = gff_df.loc[gff_df[2]=="insertion_sequence",:]

    for idx,row in sub_df.iterrows():
        IS_info = dict([_.split('=') for _ in row[8].split(';')])
        IS_ID = IS_info['ID']
        family = IS_info.get('family','UNKNOWN')
        cluster = IS_info.get('cluster','UNKNOWN')

        contig = row[0]
        seq = contig2seq[contig]
        s,e = row[3],row[4]
        sub_seq = seq[s:e]

        all_fea = [_ for _ in sub_seq.features if 'locus_tag' in _.qualifiers]
        all_locus = [fea.qualifiers['locus_tag'][0] for fea in all_fea]
        all_kos = sorted(list(set([locus2ko.get(_,[_])[0]
                                   for _ in all_locus])))
        all_ko_enzymes = ';'.join(sorted(list(set([locus2ko.get(_,[_,_])[1]
                                                   for _ in all_locus]))))
       # IS_ID2locus_list[IS_ID] = [locus2ko.get(_,_) for _ in all_locus]
        IS_info_df.loc[IS_ID,:] = [cluster,s,e,row[6],
                                   family,
                                   len(sub_seq.seq),
                                   len(all_fea),
                                   len(all_kos),
                                   all_ko_enzymes]
IS_info_df.loc[:,'genome'] = [_.split('_')[0] for _ in IS_info_df.index]
IS_info_df.loc[:,'contig'] = [_.rsplit('_',2)[0] for _ in IS_info_df.index]
IS_info_df.loc[:,'MC name'] = [g2pop[_] for _ in IS_info_df['genome']]
# IS_info_df.loc[:,'in chrom'] = ['Yes'
#                                 if _ in g2chrom.get(_.split('_')[0],[]) else 'No' for _ in IS_info_df['contig']]

sub = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/pub_dataset/IS_info.tsv',sep='\t',index_col=0)
IS_info_df = pd.concat([IS_info_df,sub],axis=0)
IS_info_df.to_csv(f"{basedir}/canu_o/IS_info.tsv",sep='\t',index=0)






    




