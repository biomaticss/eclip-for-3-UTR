import pandas as pd

# 配置
INPUT_TXT = 'human.PAS.txt'  # PolyA_DB输入文件
OUTPUT_BED = 'target_pas_hg19.bed'  # 输出BED文件
GENES = ['STAT5B', 'GSPT1', 'CCND2', 'CTNNBIP1']  # 目标基因

# 转换函数
def convert_pas_to_bed():
    # 读取human.PAS.txt
    df = pd.read_csv(INPUT_TXT, sep='\t', low_memory=False)
    
    # 筛选目标基因
    df = df[df['Gene Symbol'].isin(GENES)]
    
    # 筛选3' UTR或3' most exon的PAS
    df = df[df['Intron/exon location'].str.contains("3'|UTR", case=False, na=False)]
    
    # 创建BED格式
    bed_df = pd.DataFrame({
        'chrom': df['Chromosome'],
        'start': df['Position'] - 1,  # 0-based start
        'end': df['Position'],       # 1-based end
        'name': df['PAS_ID'],
        'score': df['Mean RPM'],
        'strand': df['Strand'],
        'gene_symbol': df['Gene Symbol']  # 附加列，便于后续筛选
    })
    
    # 保存为BED（6列+额外gene_symbol）
    bed_df.to_csv(OUTPUT_BED, sep='\t', header=False, index=False)
    print(f'Generated {OUTPUT_BED} with {len(bed_df)} PAS entries for {GENES}')

if __name__ == '__main__':
    convert_pas_to_bed()