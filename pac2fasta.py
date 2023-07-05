#!/usr/bin/env python
# coding: utf-8
import argparse
from collections import namedtuple

parser = argparse.ArgumentParser(
    prog='pac2fasta',
    description='get fasta from pac,ann,amb',
    )

parser.add_argument('fasta_path', help='fasta path for output fasta and input index files')
args = parser.parse_args()

target_fa = args.fasta_path
index_pac = target_fa+".pac"
index_ann = target_fa+".ann"
index_amb = target_fa+".amb"

ANN = namedtuple("ANN","gi name anno offset len amb_num")
ann_list = []
with open(index_ann, 'r') as f:
    f.readline()
    while True:
        line1 = f.readline()
        line2 = f.readline()
        if not line2:
            break
        _l1 = line1.strip().split()
        _l2 = line2.strip().split()
        _gi, _name, _anno, _offset, _len, _amb_num = _l1+_l2
        _offset = int(_offset)
        _len = int(_len)
        _amb_num = int(_amb_num)
        _anno = '' if _anno == '(null)' else _anno
        ann_list.append(ANN(_gi, _name, _anno, _offset, _len, _amb_num))

AMB = namedtuple("AMB","offset len base")
amb_list = []
with open(index_amb, 'r') as f:
    f.readline()
    for line in f:
        _offset, _len, _base = line.strip().split()
        _offset = int(_offset)
        _len = int(_len)
        amb_list.append(AMB(_offset, _len, _base))

with open(index_pac, 'br') as f:
    bs = f.read()
    binary_string = bin(int(bs.hex(), 16))
    bss = binary_string[2:]

base_dict = {"11": "T", "01": "C", "10": "G", "00": "A"}

ConcatSeq = [base_dict["".join(bss[i:i+2])] for i in range(0, len(bss), 2)]
    
for amb in amb_list:
    ConcatSeq[amb.offset:amb.offset+amb.len] = [amb.base] * amb.len

with open(target_fa, 'w') as f:
    for ann in ann_list:
        _header = f">{ann.name} {ann.anno}\n"
        _seq = ConcatSeq[ann.offset:ann.offset+ann.len]
        _seq_80 = ["".join(_seq[i:i+80]) for i in range(0, len(_seq), 80)]
        f.write(_header)
        f.write("\n".join(_seq_80)+'\n')

