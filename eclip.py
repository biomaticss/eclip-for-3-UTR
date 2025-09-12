import pybedtools
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns

# 尝试导入pyBigWig，如果失败则给出提示
try:
    import pyBigWig
    PYBIGWIG_AVAILABLE = True
except ImportError:
    PYBIGWIG_AVAILABLE = False
    print("Warning: pyBigWig module not found. Signal extraction will be skipped.")
    print("To install: pip install pyBigWig")

# 配置
GENES = ['STAT5B', 'GSPT1', 'CCND2', 'CTNNBIP1']
RBP_FILES = {
    'PTBP1': {
        'bed': 'PTBP1_K562.bed',
        'bigwig_plus': 'PTBP1_K562_plus.bigWig',
        'bigwig_minus': 'PTBP1_K562_minus.bigWig',
        'cell_line': 'K562'
    },
    'DDX3X': {
        'bed': 'DDX3X_HepG2.bed',
        'bigwig_plus': 'DDX3X_HepG2_plus.bigWig',
        'bigwig_minus': 'DDX3X_HepG2_minus.bigWig',
        'cell_line': 'HepG2'
    },
    'ELAVL1': {
        'bed': 'ELAVL1_K562.bed',
        'bigwig_plus': 'ELAVL1_K562_plus.bigWig',
        'bigwig_minus': 'ELAVL1_K562_minus.bigWig',
        'cell_line': 'K562'
    },
    'PABPN1': {
        'bed': 'PABPN1_HepG2.bed',
        'bigwig_plus': 'PABPN1_HepG2_plus.bigWig',
        'bigwig_minus': 'PABPN1_HepG2_minus.bigWig',
        'cell_line': 'HepG2'
    }
}
PAS_NEARBY_BED = 'pas_nearby.bed'
GENE_BED = 'gene_coords.bed'
CHROM_SIZES = 'hg38.chrom.sizes'
WINDOW = 500  # PAS附近窗口大小
BIN_SIZE = 50  # 信号分箱大小

# 交集分析
def intersect_peaks_with_pas(rbp, bed_file, pas_bed):
    peaks = pybedtools.BedTool(bed_file)
    pas_regions = pybedtools.BedTool(pas_bed)
    intersections = peaks.intersect(pas_regions, wa=True, wb=True)
    output_file = f'{rbp}_intersections.bed'
    intersections.saveas(output_file)
    return output_file

# 使用pyBigWig提取bigWig信号
def extract_signal_around_pas(bigwig_file, chrom, center_pos, window=500, bin_size=50):
    """提取PAS周围区域的信号"""
    if not PYBIGWIG_AVAILABLE:
        return [], []
    
    try:
        if not os.path.exists(bigwig_file):
            print(f"Warning: bigWig file {bigwig_file} not found")
            return [], []
        
        bw = pyBigWig.open(bigwig_file)
        if not bw.isBigWig():
            print(f"Warning: {bigwig_file} is not a valid bigWig file")
            bw.close()
            return [], []
        
        start = max(0, center_pos - window)
        end = center_pos + window
        
        # 获取信号值
        signal = bw.values(chrom, start, end)
        bw.close()
        
        if signal is None:
            return [], []
        
        # 将信号分箱
        bins = []
        bin_values = []
        for i in range(0, len(signal), bin_size):
            bin_start = start + i
            bin_end = min(end, bin_start + bin_size)
            # 处理可能的NaN值
            bin_chunk = signal[i:i+bin_size]
            # 过滤掉NaN值
            valid_values = [x for x in bin_chunk if x is not None and not np.isnan(x)]
            if valid_values:
                bin_val = np.mean(valid_values)
            else:
                bin_val = 0
            bins.append((bin_start + bin_end) // 2)
            bin_values.append(bin_val)
            
        return bins, bin_values
    except Exception as e:
        print(f"Error reading {bigwig_file}: {e}")
        return [], []

# 可视化函数
def plot_heatmap(rbps, gene_bed, pas_bed, output_dir='plots'):
    """创建热图展示RBP在PAS附近的结合信号"""
    if not PYBIGWIG_AVAILABLE:
        print("Skipping heatmap generation: pyBigWig not available")
        return
        
    os.makedirs(output_dir, exist_ok=True)
    gene_df = pd.read_csv(gene_bed, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    pas_df = pd.read_csv(pas_bed, sep='\t', header=None, 
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'gene'],
                        dtype={'gene': str})
    
    for gene in GENES:
        gene_info = gene_df[gene_df['name'].str.contains(gene, case=False, na=False)]
        if gene_info.empty:
            print(f"No gene coordinates for {gene}")
            continue
            
        gene_info = gene_info.iloc[0]
        chrom, gstrand = gene_info['chrom'], gene_info['strand']
        
        pas_info = pas_df[pas_df['gene'].str.contains(gene, case=False, na=False)]
        if pas_info.empty:
            print(f"No PAS found for {gene}")
            continue
            
        pas_pos = (pas_info['start'].iloc[0] + pas_info['end'].iloc[0]) // 2
        
        # 准备热图数据
        heatmap_data = []
        rbp_names = []
        
        for rbp, info in rbps.items():
            if gstrand == '+':
                bigwig_file = info.get('bigwig_plus')
            else:
                bigwig_file = info.get('bigwig_minus')
                
            # 检查bigwig文件是否存在
            if not bigwig_file:
                continue
                
            bins, signals = extract_signal_around_pas(bigwig_file, chrom, pas_pos, WINDOW, BIN_SIZE)
            if signals:
                heatmap_data.append(signals)
                rbp_names.append(rbp)
        
        if not heatmap_data:
            print(f"No signal data for {gene}")
            continue
            
        # 创建热图
        plt.figure(figsize=(12, 8))
        sns.heatmap(heatmap_data, xticklabels=len(bins)//5 if bins else 1, yticklabels=rbp_names, 
                   cmap='YlOrRd', cbar_kws={'label': 'Signal intensity'})
        
        plt.xlabel('Genomic Position (bp)')
        plt.ylabel('RBP')
        plt.title(f'RBP Binding Signal Heatmap around PAS for {gene} (strand: {gstrand})')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()

def plot_bar_chart(rbps, gene_bed, pas_bed, output_dir='plots'):
    """创建柱状图展示每个RBP在PAS附近的结合强度"""
    if not PYBIGWIG_AVAILABLE:
        print("Skipping bar chart generation: pyBigWig not available")
        return
        
    os.makedirs(output_dir, exist_ok=True)
    gene_df = pd.read_csv(gene_bed, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    pas_df = pd.read_csv(pas_bed, sep='\t', header=None, 
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'gene'],
                        dtype={'gene': str})
    
    for gene in GENES:
        gene_info = gene_df[gene_df['name'].str.contains(gene, case=False, na=False)]
        if gene_info.empty:
            print(f"No gene coordinates for {gene}")
            continue
            
        gene_info = gene_info.iloc[0]
        chrom, gstrand = gene_info['chrom'], gene_info['strand']
        
        pas_info = pas_df[pas_df['gene'].str.contains(gene, case=False, na=False)]
        if pas_info.empty:
            print(f"No PAS found for {gene}")
            continue
            
        pas_pos = (pas_info['start'].iloc[0] + pas_info['end'].iloc[0]) // 2
        
        # 准备柱状图数据
        rbp_names = []
        avg_signals = []
        
        for rbp, info in rbps.items():
            if gstrand == '+':
                bigwig_file = info.get('bigwig_plus')
            else:
                bigwig_file = info.get('bigwig_minus')
                
            # 检查bigwig文件是否存在
            if not bigwig_file:
                continue
                
            bins, signals = extract_signal_around_pas(bigwig_file, chrom, pas_pos, WINDOW, BIN_SIZE)
            if signals:
                avg_signal = np.mean(signals)
                rbp_names.append(rbp)
                avg_signals.append(avg_signal)
        
        if not avg_signals:
            print(f"No signal data for {gene}")
            continue
            
        # 创建柱状图
        plt.figure(figsize=(10, 6))
        colors = plt.cm.Set3(np.linspace(0, 1, len(rbp_names)))
        bars = plt.bar(rbp_names, avg_signals, color=colors)
        
        # 添加数值标签
        for bar, value in zip(bars, avg_signals):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                    f'{value:.2f}', ha='center', va='bottom')
        
        plt.xlabel('RBP')
        plt.ylabel('Average Signal Intensity')
        plt.title(f'Average RBP Binding Intensity around PAS for {gene} (strand: {gstrand})')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_barchart.png', dpi=300, bbox_inches='tight')
        plt.close()

def plot_binding_profile(rbps, gene_bed, pas_bed, output_dir='plots'):
    """创建结合谱线图展示每个RBP在PAS附近的信号分布"""
    if not PYBIGWIG_AVAILABLE:
        print("Skipping binding profile generation: pyBigWig not available")
        return
        
    os.makedirs(output_dir, exist_ok=True)
    gene_df = pd.read_csv(gene_bed, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    pas_df = pd.read_csv(pas_bed, sep='\t', header=None, 
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'gene'],
                        dtype={'gene': str})
    
    for gene in GENES:
        gene_info = gene_df[gene_df['name'].str.contains(gene, case=False, na=False)]
        if gene_info.empty:
            print(f"No gene coordinates for {gene}")
            continue
            
        gene_info = gene_info.iloc[0]
        chrom, gstrand = gene_info['chrom'], gene_info['strand']
        
        pas_info = pas_df[pas_df['gene'].str.contains(gene, case=False, na=False)]
        if pas_info.empty:
            print(f"No PAS found for {gene}")
            continue
            
        pas_pos = (pas_info['start'].iloc[0] + pas_info['end'].iloc[0]) // 2
        
        plt.figure(figsize=(12, 8))
        
        for rbp, info in rbps.items():
            if gstrand == '+':
                bigwig_file = info.get('bigwig_plus')
            else:
                bigwig_file = info.get('bigwig_minus')
                
            # 检查bigwig文件是否存在
            if not bigwig_file:
                continue
                
            bins, signals = extract_signal_around_pas(bigwig_file, chrom, pas_pos, WINDOW, BIN_SIZE)
            if signals and bins:
                plt.plot(bins, signals, label=f'{rbp} ({info["cell_line"]})', linewidth=2)
        
        plt.axvline(x=pas_pos, color='red', linestyle='--', label='PAS position')
        plt.xlabel('Genomic Position (bp)')
        plt.ylabel('Signal Intensity')
        plt.title(f'RBP Binding Profile around PAS for {gene} (strand: {gstrand})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_profile.png', dpi=300, bbox_inches='tight')
        plt.close()

# 主程序
def main():
    # 创建交集文件
    for rbp, info in RBP_FILES.items():
        print(f'Processing {rbp}...')
        if os.path.exists(info['bed']):
            intersect_peaks_with_pas(rbp, info['bed'], PAS_NEARBY_BED)
        else:
            print(f"Warning: {info['bed']} does not exist")
    
    # 创建多种可视化
    print("Creating heatmaps...")
    plot_heatmap(RBP_FILES, GENE_BED, PAS_NEARBY_BED)
    
    print("Creating bar charts...")
    plot_bar_chart(RBP_FILES, GENE_BED, PAS_NEARBY_BED)
    
    print("Creating binding profiles...")
    plot_binding_profile(RBP_FILES, GENE_BED, PAS_NEARBY_BED)
    
    if not PYBIGWIG_AVAILABLE:
        print("\nNote: Signal visualization was skipped because pyBigWig is not available.")
        print("Peak intersection analysis was still performed.")

if __name__ == '__main__':
    main()