import os
import random
import glob
import argparse
import threading
import logging
import json
from scripts import ont_qc_info

# 安静删除文件和目录
def delete_files(files):
    """
    files: [ 'file1', 'file2' , ...]
    """
    try:
        for f in files:
            try:
                os.remove(f)
            except:
                pass
    except:
        pass

def randon_code(n=6,alpha=True):
    """
    randon_code(6,False)  # 打印6位数字验证码
    randon_code(4,True)   # 打印4位数字字母混合验证码
    """
    s = '' # 创建字符串变量,存储生成的验证码
    for i in range(n):  # 通过for循环控制验证码位数
        num = random.randint(0,9)  # 生成随机数字0-9
        if alpha: # 需要字母验证码,不用传参,如果不需要字母的,关键字alpha=False
            upper_alpha = chr(random.randint(65,90))
            lower_alpha = chr(random.randint(97,122))
            num = random.choice([num,upper_alpha,lower_alpha])
        s = s + str(num)
    return s

def getFileSize(filePath):
    try:
        fsize = os.path.getsize(filePath)   # 返回的是字节大小
        """
        为了更好地显示，应该时刻保持显示一定整数形式，即单位自适应
        return { 'size': xxxx, 'size_f': 'xxxG'}
        """
        if fsize < 1024:
            return { 'size': fsize, 'size_f': str(round(fsize,2)) + 'Byte'}
        else:
            KBX = fsize/1024
            if KBX < 1024:
                return { 'size': fsize, 'size_f': str(round(KBX,2)) + 'K'}
            else:
                MBX = KBX /1024
                if MBX < 1024:
                    return { 'size': fsize, 'size_f': str(round(MBX,2)) + 'M'}
                else:
                    return { 'size': fsize, 'size_f': str(round(MBX/1024)) + 'M'}
    except:
        return { 'size': -1, 'size_f': '' }

def file_check(filePath):
    """
    input : a file path
    return : { 'file': filePath, 'size': size, 'isExist': True}
    """
    filePath = filePath.replace('\>', '>')
    return { 'file': filePath, 'size': getFileSize(filePath)['size'], 'size_f': getFileSize(filePath)['size_f'], 'isExist': os.path.exists(filePath)}

def find_fqgz(absPath):
    """
    input: 
    absPath: /pfs/byProject/WHXWZM-2022120162A/SequenceData/Nanopore/1D/TGXY2300126802/20231023-NPL2300376-P4-PAS48514-sup/PAS48514/
    return: { 'status':'成功', 'data_pass': [pass.fq.gz1, pass.fq.gz2, ...],'data_fail': [fail.fq.gz1, fail.fq.gz2, ...], 'info':err }
    """
    try:
        # 找qc_report下的fq.gz和summary
        absPath = absPath.replace('\>', '>')
        data_pass = glob.glob(absPath + '/*qc_report/*pass*fastq.gz') + glob.glob(absPath + '/*qc_report/barcode*/*pass*fastq.gz')
        data_fail = glob.glob(absPath + '/*qc_report/*fail*fastq.gz') + glob.glob(absPath + '/*qc_report/barcode*/*fail*fastq.gz')
        summary = glob.glob(absPath + '/sequencing_summary/*summary*txt')
        return { 'status':'成功', 'data_pass': data_pass,'data_fail': data_fail, 'info':'', 'summary': summary }
    except Exception as e:
        return { 'status':'失败', 'data_pass': [],'data_fail': [], 'info':e, 'summary':[] }

def find_bam(absPath):
    """
    input: 
    absPath: /pfs/byProject/WHXWZM-2022120162A/SequenceData/Nanopore/1D/TGXY2300126802/20231023-NPL2300376-P4-PAS48514-sup/PAS48514/
    return: { 'status':'成功', 'data_pass': [pass.fq.gz1, pass.fq.gz2, ...],'data_fail': [fail.fq.gz1, fail.fq.gz2, ...], 'info':err }
    """
    try:
        # 找qc_report下的bam和summary
        data_pass = glob.glob(absPath + '/bam_pass/*pass*bam') + glob.glob(absPath + '/bam_pass/barcode*/*pass*bam')
        data_fail = glob.glob(absPath + '/bam_fail/*fail*bam') + glob.glob(absPath + '/bam_fail/barcode*/*fail*bam')
        summary = glob.glob(absPath + '/bam_fail/*fail.summary.txt') + glob.glob(absPath + '/bam_pass/*pass.summary.txt')
        return { 'status':'成功', 'data_pass': data_pass,'data_fail': data_fail, 'info':'' , 'summary': summary }
    except Exception as e:
        return { 'status':'失败', 'data_pass': [],'data_fail': [], 'info':e, 'summary': [] }

# 合并ONT fastq.gz/bam
def merge_ont(file_list, outfile):
    """
    input:
    data_list: [ file1, file2, file3 ...]
    outfile: merge file1 file2 file3 ... to outfile
    return: 
    { 'status':'完成', 'data': 'outfile abs path', 'info':'some logs', 'size': '??G' }
    """
    new_source_paths = []
    for i in file_list:
        i = i.replace('>', '\>')
        new_source_paths.append(i)
    file_list = new_source_paths
    # print(file_list)
    try:
        assert len(file_list)
        # 检查提供的的文件是否都存在
        for filePath in file_list:
            file = file_check(filePath)
            if not file['isExist']:
                raise Exception(f'文件 {filePath} 不存在')
        # 创建临时合并目录和合并
        work_dir  = os.path.dirname(outfile)
        file_name = os.path.basename(outfile)
        dtype = 'pass' if 'pass' in file_name else 'fail'
        if 'fast' in os.path.basename(file_list[0]):
            ddtype = 'fast.' + dtype
        elif 'sup' in os.path.basename(file_list[0]):
            ddtype = 'sup.' + dtype
        elif 'hac' in os.path.basename(file_list[0]):
            ddtype = 'hac.' + dtype
        else:
            ddtype = dtype
        file_name = file_name.replace(dtype, ddtype)
        outfile = os.path.join(work_dir, file_name)
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        try:
            if file_list[0].endswith('bam'):
                # 开始合并
                print(f'开始合并{dtype}.fastq.gz')
                outfile_temp1 = os.path.join(work_dir, str(file_name) + '.' + randon_code(6,True))
                os.system(f"samtools merge --threads 8  -o {outfile_temp1} {' '.join(file_list)}  && mv {outfile_temp1} {outfile}")
                res = file_check(outfile)
            elif file_list[0].endswith('fastq.gz'):
                # 开始合并
                print(f'开始合并{dtype}.fastq.gz')
                outfile_temp1 = os.path.join(work_dir, str(file_name) + '.' + randon_code(6,True))
                os.system(f"cat {' '.join(file_list)} > {outfile_temp1} && mv {outfile_temp1} {outfile}")
                res = file_check(outfile)
            else:
                raise Exception('不符合数据格式 fastq.gz/bam')
        except:
            raise Exception('合并失败！')
        return { 'status':'完成', 'data': outfile, 'info':'', 'size': res['size'], 'size_f': res['size_f']  }
    except Exception as e:
        err = '|'.join([e.__traceback__.tb_frame.f_globals["__file__"], str(e.__traceback__.tb_lineno), str(e)])
        print(json.dumps({ 'status':'失败', 'data': '', 'info':err, 'size': '', 'size_f': '' }))
        return { 'status':'失败', 'data': '', 'info':err, 'size': '', 'size_f': '' }

def qc_combined_data(pass_fq, fail_fq, summary):
    return 0

class args1:
    fq_lst = ''
    summary_lst = ''
    out_dir = ''
    def __init__(self,fq_lst,summary_lst,out_dir):
        self.fq_lst = fq_lst
        self.summary_lst = summary_lst
        self.out_dir = out_dir

# ONT 质控
def ont_qc(files, summarys, workdir):
    logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(asctime)s %(message)s')
    LOG = logging.getLogger(__name__)
    info_path = os.path.join(workdir, 'ont_qc_info')
    args = args1(files, summarys, info_path)
    LOG.info(args)
    ont_qc_info.normal_run(args, LOG)   # ont qc info
    info_list = os.path.join(workdir, 'info.lst')
    os.system(f'realpath {info_path}/*txt > {info_list}')
    ont_qc_info.StatAndFig.create_stats_and_figs(info_list, workdir)  # json + 3 pngs

## 合并 ONT 数据
def ont_deal_input(source_paths, out_path, name_prefix, qc):
    # 处理输入路径
    new_source_paths = []
    for i in source_paths:
        i = i.replace('>', '\>')
        new_source_paths.append(i)
    source_paths = new_source_paths
    # raise Exception('')
    # 存放找到的待合并文件
    datas = dict()
    datas['fq_pass'] = []
    datas['fq_fail'] = []
    datas['bam_pass'] = []
    datas['bam_fail'] = []
    datas['summary'] = []
    # 输出文件名
    fastq_pass_file = os.path.join(out_path, name_prefix + '.pass.fastq.gz')
    fastq_fail_file = os.path.join(out_path, name_prefix + '.fail.fastq.gz')
    bam_pass_file = os.path.join(out_path, name_prefix + '.pass.bam')
    bam_fail_file = os.path.join(out_path, name_prefix + '.fail.bam')
    # 根据路经查找数据
    for path in source_paths:
        datas['fq_pass'] += find_fqgz(path)["data_pass"]
        datas['fq_fail'] += find_fqgz(path)["data_fail"]
        datas['bam_pass'] += find_bam(path)["data_pass"]
        datas['bam_fail'] += find_bam(path)["data_fail"]
        datas['summary'] += find_fqgz(path)["summary"] + find_bam(path)["summary"]
    # 合并
    threads = []
    print(datas)
    if len(datas['fq_pass']) >= 1:
        threads.append(threading.Thread(target=merge_ont, args=(datas['fq_pass'], fastq_pass_file)))
    if len(datas['fq_fail']) >= 1:
        threads.append(threading.Thread(target=merge_ont, args=(datas['fq_fail'], fastq_fail_file)))
    if len(datas['bam_pass']) >= 1:
        threads.append(threading.Thread(target=merge_ont, args=(datas['bam_pass'], bam_pass_file)))
    if len(datas['bam_fail']) >= 1:
        threads.append(threading.Thread(target=merge_ont, args=(datas['bam_fail'], bam_fail_file)))
    # 质控
    if len(threads) > 0 and qc == 'yes':
        # 质控目录
        qc_path = os.path.join(out_path, 'ont_qc_report')
        if not os.path.exists(qc_path):
            os.mkdir(qc_path)
        summary_lst = os.path.join(qc_path, 'summary.lst')
        fq_list = os.path.join(qc_path, 'fq.list')
        open(fq_list, 'w').write('\n'.join(datas['fq_pass'] + datas['fq_fail'] + datas['bam_pass'] + datas['bam_fail']))
        open(summary_lst, 'w').write('\n'.join(datas['summary']))
        threads.append(threading.Thread(target=ont_qc, args=(fq_list, summary_lst, qc_path)))
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return 1


def main():
    parser = argparse.ArgumentParser(description="合并数据")
    parser.add_argument('--sources', '-s', type=str, nargs='+', required=True, help="输入一个或多个标准目录路经, 例如: \'/pfs/byProject/WHXWZM-2022120162A/SequenceData/Nanopore/1D/TGXY2300126802/20231023-NPL2300376-P4-PAS48514-sup/PAS4851\'")
    parser.add_argument('--out', '-o', type=str, required=True, help="输出合并后的文件存放目录")
    parser.add_argument('--prefix', '-p', type=str, required=True, help="输出合并后的文件名前缀")
    parser.add_argument('--qc', '-q', type=str, choices=['yes', 'no'], default='no', help="是否质控(可能会质控失败，但不影响合并) 默认:no")
    parser.add_argument('--log', '-l', type=str, choices=['yes', 'no'], default='no', help="是否输出合并的日志文件 默认:no")
    args = parser.parse_args()
    ont_deal_input(args.sources, args.out, args.prefix, args.qc)

if __name__ == "__main__":
    main()





