import pandas as pd
import subprocess
import os

# 配置
GENES = ['STAT5B', 'GSPT1', 'CCND2', 'CTNNBIP1']
HG38_PAS_BED = 'target_pas_hg38.bed'  # liftOver后的hg38 PAS BED
PROXIMAL_PAS_BED = 'proximal_pas.bed'  # 近端PAS输出
PAS_NEARBY_BED = 'pas_nearby.bed'  # 扩展±200 bp
CHROM_SIZES = 'hg38.chrom.sizes'  # 从UCSC下载

# 终止密码子位置（hg38，需确认）
STOP_CODON = {
    'STAT5B': 42276325,  # +链, CDS end
    'GSPT1': 11877048,   # -链, CDS start
    'CCND2': 4302728,    # +链, CDS end
    'CTNNBIP1': 11026219 # -链, CDS start
}

# 筛选近端PAS
def select_proximal_pas():
    # 读取PAS BED
    df = pd.read_csv(HG38_PAS_BED, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'pas_id', 'score', 'strand', 'gene'])
    
    proximal = []
    for gene in GENES:
        gene_df = df[df['gene'] == gene]
        if gene_df.empty:
            print(f"No PAS found for {gene}")
            continue
        
        stop = STOP_CODON[gene]
        gene_strand = gene_df['strand'].iloc[0]
        
        if gene_strand == '+':
            # 正链：PAS end > stop codon, 取最小的end
            gene_df = gene_df[gene_df['end'] > stop]
            if gene_df.empty:
                print(f"No valid 3' UTR PAS for {gene}")
                continue
            prox = gene_df.loc[gene_df['end'].idxmin()]
        else:
            # 负链：PAS start < stop codon, 取最大的start
            gene_df = gene_df[gene_df['start'] < stop]
            if gene_df.empty:
                print(f"No valid 3' UTR PAS for {gene}")
                continue
            prox = gene_df.loc[gene_df['start'].idxmax()]
        
        proximal.append(prox)
    
    # 保存近端PAS
    if proximal:
        prox_df = pd.DataFrame(proximal)
        prox_df.to_csv(PROXIMAL_PAS_BED, sep='\t', header=False, index=False, 
                       columns=['chrom', 'start', 'end', 'pas_id', 'score', 'strand'])
        print(f'Generated {PROXIMAL_PAS_BED} with {len(prox_df)} proximal PAS')
    else:
        print("No proximal PAS found")
        return

    # 扩展±200 bp
    cmd = f'bedtools slop -i {PROXIMAL_PAS_BED} -g {CHROM_SIZES} -b 200 > {PAS_NEARBY_BED}'
    subprocess.run(cmd, shell=True, check=True)
    print(f'Generated {PAS_NEARBY_BED}')

if __name__ == '__main__':
    select_proximal_pas()