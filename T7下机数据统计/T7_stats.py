# -*- coding: utf-8 -*-
# @Author  : Gongshuang
# @Time    : 2023/4/10 11:00
# @Function: Stat T7 data
# @Document: https://grandomics.feishu.cn/docx/FXwPd8yf9o8NBlx2bChccGcnnrb?from=from_copylink
import argparse
import sys
import os
import logging
import time
import datetime
import subprocess

LOG = logging.getLogger()

# 下机数据根路径
ROOT = '/pfs/Sequencing/DNBSEQ-T7/rawData/R1100600210024'


def check_path(path):
    if not os.path.exists(path):
        LOG.error("No such file or directory: {}".format(path))
        raise Exception

def get_statfq_info(fin):
    dict_r = {}
    for i in fin:
        name = i.strip().split('\t')[0]
        value = i.strip().split('\t')[1]
        dict_r[name[1:]] = value
        if name == '#EstErr%':
           break
    return dict_r

def count_total(r1, r2, b1, b2):
    r1 = float(r1)/100.0
    r2 = float(r2)/100.0
    b1 = int(b1)
    b2 = int(b2)
    return "%.2f" % ((r1*b1 + r2*b2)/(b1 + b2) * 100)

# 解析R1 和 R2 的fqStat.txt文件，得到：reads bases(G) GC(%) Q20(%) Q30(%) 补测/G
def parse_fqStat(fqStat1, fqStat2):
    dic_r1 = {}
    dic_r2 = {}
    try:
        with open(fqStat1) as fin1, open(fqStat2) as fin2:
            dic_r1 = get_statfq_info(fin1)
            dic_r2 = get_statfq_info(fin2)
        total_GC = count_total(dic_r1['GC%'], dic_r2['GC%'], dic_r1['BaseNum'], dic_r2['BaseNum'])
        total_Q20 = count_total(dic_r1['Q20%'], dic_r2['Q20%'], dic_r1['BaseNum'], dic_r2['BaseNum'])
        total_Q30 = count_total(dic_r1['Q30%'], dic_r2['Q30%'], dic_r1['BaseNum'], dic_r2['BaseNum'])
        return [dic_r1['ReadNum'], int(dic_r1['BaseNum'])*2, total_GC, total_Q20, total_Q30]
    except:
        return [0, 0, 0, 0, 0]


def resolve_files(source_dir, sno, barcodes, require, outfile):
    source_dir = os.path.abspath(source_dir)
    listdir = os.listdir(source_dir)
    # 文件前缀
    prefix = [s[:14] for s in listdir if 'fq.fqStat.txt' in s][0]
    total_bases = 0
    barcode_stat = {}  # 每个barcode有自己的统计行
    for barcode in barcodes:
        Read1 = source_dir + '/{prefix}_{barcode}_1.fq.gz'.format(prefix=prefix, barcode=barcode)
        Read2 = source_dir + '/{prefix}_{barcode}_2.fq.gz'.format(prefix=prefix, barcode=barcode)
        fqStat1 = source_dir + '/{prefix}_{barcode}_1.fq.fqStat.txt'.format(prefix=prefix, barcode=barcode)
        fqStat2 = source_dir + '/{prefix}_{barcode}_2.fq.fqStat.txt'.format(prefix=prefix, barcode=barcode)
        reads, bases, GC, Q20, Q30 = parse_fqStat(fqStat1, fqStat2)
        total_bases += int(bases)  # 计算一个样本多个barcode的碱基总数，用于计算是否需要补测
        barcode_stat[barcode] = [reads, bases, GC, Q20, Q30, Read1, Read2]
        LOG.info("Counting for barcode : {}".format(barcode))
    # barcode_stat[barcodes[0]][1] = total_bases
    #统计补测数据
    buce = "%.2f" % (float(require) - total_bases/1000000000.0)
    if float(buce) <= 0:
        buce = 0
    barcode_stat[barcodes[0]].append(buce)
    # 统计结果写入文件
    # LOG.info("Writing to file : T7_stat_Result.xls")
    for barcode in barcodes:
        reads = int(barcode_stat[barcode][0])
        bases = '%.2f' % (float(barcode_stat[barcode][1])/1000000000.0)
        GC = barcode_stat[barcode][2]
        Q20 = barcode_stat[barcode][3]
        Q30 = barcode_stat[barcode][4]
        Read1 = barcode_stat[barcode][5]
        Read2 = barcode_stat[barcode][6]
        try:
            buce = barcode_stat[barcode][7]
        except:
            buce = ''
        outfile.write("{}\t{}\t{}\t{:,}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sno, barcode,require ,reads, str(bases), GC, Q20, Q30, str(buce), Read1, Read2))

def T7_stats(source_dir, input, out_dir, outfile):
    check_path(source_dir)
    check_path(out_dir)
    input = os.path.abspath(input)
    check_path(input)
    LOG.info("Reading directory: {}".format(source_dir))
    LOG.info("Reading file: {}".format(input))
    # reading barcode file
    bcInfo = []
    bcInfoDict = {}
    with open(input) as bc:
        for i in bc:
            sno, barcode, require = i.strip().split('\t')
            bcInfo.append([sno, barcode, require])
            bcInfoDict[sno] = []
        # print(bcInfoDict)
        for line in bcInfo:
            bcInfoDict[line[0]].append(line[1])
        #print(bcInfoDict)
    counted_sno = []
    for line in bcInfo:
        if line[0] in counted_sno:
            continue
        #do_sth(line[0], bcInfoDict[line[0]], line[2])
        resolve_files(source_dir, line[0], bcInfoDict[line[0]], line[2],  outfile)
        counted_sno.append(line[0])

# 返回shell命令结果
def file_process(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8', executable='/bin/bash')
    out, e = p.communicate()
    return_code = p.returncode
    if return_code == 0:
        if out != '':
            return out
        return True
    else:
        return False

def main():
    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                        format="[%(levelname)s] %(asctime)s %(message)s")
    parser = argparse.ArgumentParser(description='To Prepare T7 stat paths')
    parser.add_argument('-c','--cell', type=str,nargs='?', default='', help='Cell id, example: E100065731')
    parser.add_argument('-i','--input', type=str, help='Barcode list file, with 3 col, sample_NO, barcode, and seq_require')
    parser.add_argument('-o','--outpath', type=str, default='./', help='Output directory')
    parser.add_argument('-p','--datapath', type=str, nargs='?', default='', help='specify only one data dir')
    args = parser.parse_args()
    if args.datapath != '' and not os.path.exists(args.datapath):
        LOG.info(f"No such file or directory : {args.datapath}")
        return 
    # 当前时间
    now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    # 输出文件
    if args.cell != '':
        out_file = os.path.join(args.outpath + 'Stat_{}-{}.xls'.format(args.cell, now_time))
    else:
        out_file = os.path.join(args.outpath + 'Stat_{}.xls'.format(now_time))
    # 根据芯片号寻找数据路径
    if args.cell != '':
        paths = file_process('realpath {ROOT}/{cell}/*/{cell}_*/'.format(ROOT=ROOT,cell=args.cell))
        paths = paths.strip().split('\n') 
        if len(paths) > 1:
            LOG.info("[WARNNING] Found more than 1 paths: {} .Please use '--datapath' or '-p' to specify only one dir.".format(paths, out_file))
            return 0
        else:
            source_dir = paths[0]
    if args.datapath != '' and os.path.exists(args.datapath):
        source_dir = args.datapath
    with open(out_file, 'w') as outfile:
        outfile.write('#样本编号\tbarcode号\t待测量\treads\tbases(G)\tGC(%)\tQ20(%)\tQ30(%)\t补测/G\t数据路径\n')
        T7_stats(source_dir, args.input, args.outpath, outfile)
    LOG.info("All barcodes complete! Results in : {}".format(out_file))
 
if __name__ == '__main__':
    main()
