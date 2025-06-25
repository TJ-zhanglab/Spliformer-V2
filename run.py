import argparse
import pdb
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from pyfaidx import Fasta
import logging
import pysam
from model import SegmentNT
from transformers import AutoTokenizer
import os
# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# os.environ["CUDA_VISIBLE_DEVICES"] =  "1"
# from multiprocessing import cpu_count

# cpu_num = 8#cpu_count() # 自动获取最大核心数目
# os.environ ['OMP_NUM_THREADS'] = str(cpu_num)
# os.environ ['OPENBLAS_NUM_THREADS'] = str(cpu_num)
# os.environ ['MKL_NUM_THREADS'] = str(cpu_num)
# os.environ ['VECLIB_MAXIMUM_THREADS'] = str(cpu_num)
# os.environ ['NUMEXPR_NUM_THREADS'] = str(cpu_num)
# torch.set_num_threads(cpu_num)


class Annotator:
    def __init__(self, ref_fasta, annotations):

        df = pd.read_csv(annotations, sep='\t', dtype={'CHROM': object})
        self.genes = df['#NAME'].to_numpy()
        self.chroms = df['CHROM'].to_numpy()
        self.strands = df['STRAND'].to_numpy()
        self.tx_starts = df['TX_START'].to_numpy()
        self.tx_ends = df['TX_END'].to_numpy()
        self.exon_starts = [np.asarray([int(i) for i in c.split(',') if i])
                            for c in df['EXON_START'].to_numpy()]
        self.exon_ends = [np.asarray([int(i) for i in c.split(',') if i])
                          for c in df['EXON_END'].to_numpy()]

        self.ref_fasta = Fasta(ref_fasta, rebuild=False)

    def get_name_and_strand(self, chrom, pos):

        chrom = normalise_chrom(chrom, list(self.chroms)[0])
        idxs = np.intersect1d(np.nonzero(self.chroms == chrom)[0],
                              np.intersect1d(np.nonzero(self.tx_starts <= pos)[0],
                                             np.nonzero(pos <= self.tx_ends)[0]))

        if len(idxs) >= 1:
            return self.genes[idxs], self.strands[idxs], idxs
        else:
            return [], [], []

    def get_pos_data(self, idx, pos):

        dist_tx_start = self.tx_starts[idx] - pos
        dist_tx_end = self.tx_ends[idx] - pos
        dist_exon_bdry = min(np.union1d(self.exon_starts[idx], self.exon_ends[idx]) - pos, key=abs)
        dist_ann = (dist_tx_start, dist_tx_end, dist_exon_bdry)
        return dist_ann

##maybe delete-2025-6-11
# def one_hot(seqs, strand):
#     # {'N':0, A': 1, 'C': 2, 'G': 3, 'T': 4}
#     def preprocess_inputs(seq):
#         fun = lambda seq: [token2int[x] for x in seq]
#         return np.transpose(np.array([[fun(seq)]]), (0, 2, 1))

#     # 将y轴1与z轴2进行交换
#     if strand == '+':
#         token2int = {x: i for i, x in enumerate('NACGT')}
#         seqs = preprocess_inputs(seqs)
#     elif strand == '-':
#         token2int = {x: i for i, x in enumerate('NTGCA')}
#         seqs = preprocess_inputs(seqs[::-1])  # 取反后互换
#     # 同样对标签进行转换，与数据集相对应
#     return seqs
def reverse_complement(seq):
    if pd.isna(seq):
        return seq
    complement = str.maketrans('ATCG', 'TAGC')
    return seq[::-1].translate(complement)

def normalise_chrom(source, target):
    def has_prefix(x):
        return x.startswith('chr')

    if has_prefix(source) and not has_prefix(target):
        return source.strip('chr')
    elif not has_prefix(source) and has_prefix(target):
        return 'chr' + source

    return source


def get_delta_scores(record, ann, mask, model, tokenizer, genotype, device):
    wid = 720
    delta_scores = []
    
    try:
        record.chrom, record.pos, record.ref, len(record.alts)
    except TypeError:
        logging.warning('Please check your variant input: {}'.format(record))
        delta_scores.append('"Input error, check your input info"')
        return delta_scores

    (genes, strands, idxs) = ann.get_name_and_strand(record.chrom, record.pos)
    if len(idxs) == 0:
        logging.warning('this position is not included in gene annotation transcript: {}'.format(record))
        delta_scores.append('"This position is not included in gene annotation transcript"')
        return delta_scores

    chrom = normalise_chrom(record.chrom, list(ann.ref_fasta.keys())[0])
    try:
        seq = ann.ref_fasta[chrom][record.pos - wid // 2-1:record.pos + wid // 2-1].seq  # get seqs,range pos-401:pos+400,the first nucleictide won't extract
    except (IndexError, ValueError):
        logging.warning('Cannot be extracted from fastafile: {}'.format(record))
        delta_scores.append('"Cannot be extracted from fastafile"')
        return delta_scores

    if seq[wid // 2:wid // 2 + len(record.ref)].upper() != record.ref:  # ensure the ref base is right
        logging.warning('The REF is not same to reference genome: {}'.format(record))
        delta_scores.append('"The REF is not same to reference genome"')
        return delta_scores

    if len(seq) != wid:
        logging.warning('The variant is near chromosome end): {}'.format(record))
        delta_scores.append('"The variant is near chromosome end"')
        return delta_scores
    
    if len(record.ref) > 1:
        logging.warning('The variant is not a SNV: {}'.format(record))  
        delta_scores.append('"Only support SNVs as input"')  
        return delta_scores
    
    if 'N' in seq.upper():
        delta_scores.append('"The sequence contains unknown nucleotide, skipping"')  # the model doesn't support seqs containing 'N'
        return delta_scores
    
    for j in range(len(record.alts)):
        for i in range(len(idxs)):

            if '.' in record.alts[j] or '-' in record.alts[j] or '*' in record.alts[j]:
                continue

            if '<' in record.alts[j] or '>' in record.alts[j]:
                continue

            if len(record.ref) > 1 and len(record.alts[j]) > 1:
                continue

            if len(record.alts[j]) > 1:
                continue
            
            
            dist_ann = ann.get_pos_data(idxs[i], record.pos)
            ref_len = len(record.ref)

            x_ref = seq.upper()
            x_alt = x_ref[:wid // 2] + str(record.alts[j]) + x_ref[wid // 2 + ref_len:]
            if strands[i] == '-':
                x_ref=reverse_complement(x_ref)
                x_alt=reverse_complement(x_alt)
            if genotype:
                # pdb.set_trace()
                if record.samples[0]['GT']==(1, 0) or record.samples[0]['GT']==(0, 1):
                    x_ref_concat=x_ref+'<eos>'+x_ref
                    x_alt_concat=x_alt+'<eos>'+x_ref
                elif record.samples[0]['GT']==(1, 1) :
                    x_ref_concat=x_ref+'<eos>'+x_ref
                    x_alt_concat=x_alt+'<eos>'+x_alt
                else:
                    delta_scores.append('"Please provide the right genotype:0/1, 1/0, 1/1"')   
            else: #if no genotype provided, treat snv as homozygosis
                x_ref_concat=x_ref+'<eos>'+x_ref
                x_alt_concat=x_alt+'<eos>'+x_alt  


            x_ref = tokenizer.batch_encode_plus([x_ref_concat], return_tensors="pt", padding="max_length", max_length = 122)["input_ids"].to(device)
            x_alt = tokenizer.batch_encode_plus([x_alt_concat], return_tensors="pt", padding="max_length", max_length = 122)["input_ids"].to(device)
            
            with torch.no_grad():
                attention_mask_ref = (x_ref != tokenizer.pad_token_id).to(device)
                attention_mask_alt = (x_alt != tokenizer.pad_token_id).to(device)

                y_ref = model(x_ref,attention_mask=attention_mask_ref).logits.cpu().numpy()[:,:,:2,0]
                y_alt = model(x_alt,attention_mask=attention_mask_alt).logits.cpu().numpy()[:,:,:2,0]
            # pdb.set_trace()
            if strands[i] == '-':  # 这边调试的时候看一下维度
                y_ref = y_ref[:, ::-1]
                y_alt = y_alt[:, ::-1]
            y_ref=y_ref[:,110:-109]
            y_alt=y_alt[:,110:-109]
            y = np.concatenate([y_ref, y_alt])  # 取要的范围
            idx_pa = (y[1, :, 0] - y[0, :, 0]).argmax()  # alt.acc-ref.acc
            idx_na = (y[0, :, 0] - y[1, :, 0]).argmax()  # ref.acc - alt.acc
            idx_pd = (y[1, :, 1] - y[0, :, 1]).argmax()  # alt.don - ref.don
            idx_nd = (y[0, :, 1] - y[1, :, 1]).argmax()  # ref.don- alt.don
            mask_pa = np.logical_and((idx_pa - 250 == dist_ann[2]), mask)
            mask_na = np.logical_and((idx_na - 250 != dist_ann[2]), mask)
            mask_pd = np.logical_and((idx_pd - 250 == dist_ann[2]), mask)
            mask_nd = np.logical_and((idx_nd - 250 != dist_ann[2]), mask)
            delta_scores.append("{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}".format(
                record.alts[j],
                genes[i],
                (y[1, idx_pa, 0] - y[0, idx_pa, 0]) * (1 - mask_pa),
                (y[0, idx_na, 0] - y[1, idx_na, 0]) * (1 - mask_na),
                (y[1, idx_pd, 1] - y[0, idx_pd, 1]) * (1 - mask_pd),
                (y[0, idx_nd, 1] - y[1, idx_nd, 1]) * (1 - mask_nd),
                idx_pa - 250,
                idx_na -250,
                idx_pd -250,
                idx_nd -250))

    return delta_scores


def get_options():
    parser = argparse.ArgumentParser(description='Version: 1.3.1')
    parser.add_argument('-I', metavar='input', nargs='?', default='./input_vcf/input.vcf',
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', metavar='output', nargs='?', default='./output_vcf/output.vcf',
                        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', metavar='reference', default='./reference/hg19.fa',
                        help='path to the reference genome fasta file')
    parser.add_argument('-A', metavar='annotation', default='./reference/grch38.txt',
                        help='"grch38" (GENCODE V24lift37 canonical annotation file in '
                             'package), "grch38" (GENCODE V24 canonical annotation file in '
                             'package), or path to a similar custom gene annotation file')
    parser.add_argument('-M', metavar='mask', nargs='?', default=0,
                        type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and '
                             'unannotated acceptor/donor loss, defaults to 0')
    parser.add_argument('-T', metavar='tissue', nargs='?', default='Amygdala',
                        help='Choose tissue to predict')

    parser.add_argument('-G', metavar='genotype', nargs='?', default=0, type=int, choices=[0, 1],
                        help='Provide genotype or not, defaults to 0 (not provided)')      
    args = parser.parse_args()

    return args


def main():
    args = get_options()

    if None in [args.I, args.O, args.M]:
        logging.error('Usage: spliformer [-h] [-I [input]] [-O [output]] -R reference -A annotation '
                      '[-M [mask]] [-T [Tissue]]')
        exit()

    try:
        vcf = pysam.VariantFile(args.I)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    header = vcf.header
    header.add_line('##INFO=<ID=Spliformer_V2,Number=.,Type=String,Description="Spliformer2.0 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor usage gain (AG), acceptor usage loss (AL), donor usage gain (DG), and donor usage loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')

    try:
        output = pysam.VariantFile(args.O, mode='w', header=header)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    Modelpath='./weight2/'+str(args.T)
    tokenizer = AutoTokenizer.from_pretrained(Modelpath, trust_remote_code=True)
    model=SegmentNT.from_pretrained(Modelpath,features=['acceptor', 'donor'])
    model.to(device)   
    model.eval()

    ann = Annotator(args.R, args.A)
    for record in vcf:
        scores = get_delta_scores(record, ann, args.M, model,tokenizer, args.G, device)
        if len(scores) > 0:
            record.info['Spliformer_V2'] = scores
        output.write(record)

    vcf.close()
    output.close()


if __name__ == '__main__':
    main()
