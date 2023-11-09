# encoding=utf-8
import gzip
import sys
import io
import os
from enum import Enum
import argparse
import functools
import numpy as np
import re
import pysam
import logging
import json
import gc
import array
from collections import defaultdict
#import psutil
import copy
import traceback


# change from str to byte, for python 2 or python 3
def to_byte1(in_str):
	return in_str


def to_byte2(in_str):
	return in_str.encode('utf-8')


to_byte = to_byte1 if sys.version_info.major == 2 else to_byte2


# change from byte to str, for python 2 or python 3
def to_str1(in_str):
	return in_str


def to_str2(in_str):
	return in_str.decode('utf-8')


to_str = to_str1 if sys.version_info.major == 2 else to_str2


# for python3: change param from str to bytes
def py3_wrap_param_to_bytes(f):
	@functools.wraps(f)
	def f2(*largs, **kwargs):
		args2 = [x.encode('utf-8') if isinstance(x, str) else x for x in largs]
		return f(*args2, **kwargs)
	return f2


# for python3: change return from bytes to str
def py3_wrap_ret_to_str(f):
	@functools.wraps(f)
	def f2(*largs, **kwargs):
		ret = f(*largs, **kwargs)
		if isinstance(ret, bytes):
			return ret.decode('utf-8')
		elif isinstance(ret, (tuple, list)):
			ret2 = [x.decode('utf-8') if isinstance(x, bytes) else x for x in ret]
			return type(ret)(ret2)
		return ret
	return f2


class RWFileType(Enum):
	STD = 0x1
	GZ = 0x2
	NORMAL = 0x3


class read_file(object):

	@staticmethod
	def get_filetype_from_path(file_path):
		if file_path == '-':
			return RWFileType.STD
		elif file_path.endswith('.gz'):
			return RWFileType.GZ
		else:
			return RWFileType.NORMAL

	def __init__(self, file_path, file_type='auto', mode='r'):
		self.file_path = file_path
		self.mode = mode
		self.f = None
		if file_type == 'auto':
			self.file_type = read_file.get_filetype_from_path(file_path)
		elif file_type == '-':
			self.file_type = RWFileType.STD
		elif file_type in {'gz', 'gzip'}:
			self.file_type = RWFileType.GZ
		else:
			self.file_type = RWFileType.NORMAL

	def __enter__(self):
		if self.file_type == RWFileType.STD:
			self.f = sys.stdin
			if sys.version_info.major == 3:
				self.f.readline = py3_wrap_ret_to_str(self.f.readline)
				self.f.__iter__ = py3_wrap_ret_to_str(self.f.__iter__)
		elif self.file_type == RWFileType.GZ:
			self.f = gzip.open(self.file_path, 'rb')
			if sys.version_info.major == 3:
				self.f.readline = py3_wrap_ret_to_str(self.f.readline)
				self.f.__iter__ = py3_wrap_ret_to_str(self.f.__iter__)
		else:
			if self.mode == 'rb':
				self.f = io.open(self.file_path, mode=self.mode)
			else:
				self.f = io.open(self.file_path, mode=self.mode, encoding='utf-8')
			if sys.version_info.major == 3 and self.mode == 'rb':
				self.f.read = py3_wrap_ret_to_str(self.f.read)
		return self.f

	def __exit__(self, exc_type, exc_val, exc_tb):
		try:
			self.f.close()
		except IOError:
			pass


class write_file(object):

	def __init__(self, file_path, file_type='auto', mode='w'):
		self.file_path = file_path
		self.mode = mode
		self.f = None
		if file_type == 'auto':
			self.file_type = read_file.get_filetype_from_path(file_path)
		elif file_type == '-':
			self.file_type = RWFileType.STD
		elif file_type in {'gz', 'gzip'}:
			self.file_type = RWFileType.GZ
		else:
			self.file_type = RWFileType.NORMAL

	def __enter__(self):
		if self.file_type == RWFileType.STD:
			self.f = sys.stdout
		elif self.file_type == RWFileType.GZ:
			self.f = gzip.open(self.file_path, 'wb')
			if sys.version_info.major == 3:
				self.f.write = py3_wrap_param_to_bytes(self.f.write)
		else:
			if self.mode == 'wb':
				self.f = io.open(self.file_path, mode=self.mode)
			else:
				self.f = io.open(self.file_path, mode=self.mode, encoding='utf-8')
		return self.f

	def __exit__(self, exc_type, exc_val, exc_tb):
		try:
			self.f.close()
		except IOError:
			pass


class QC_IO(object):

	@staticmethod
	def read_lst_file(in_list_fp):
		paths = list()
		with read_file(in_list_fp) as f:
			for line in f:
				paths.append(line.strip())
		return paths

	@staticmethod
	def _read_fq(path):
		rid = seq = info = None
		with read_file(path) as f:
			ignore = False
			for line in f:
				if ignore:
					qual = line.strip()
					yield rid, seq, qual, info
					ignore = False
					continue
				if line.startswith('@'):
					info = line.strip()[1:]
					arr1 = info.split(None, 1)
					if len(arr1) >= 2:
						rid, info = arr1[0], arr1[1]
					else:
						rid, info = arr1[0], None
				elif line.startswith('+'):
					ignore = True
					continue
				else:
					seq = line.strip()

	@staticmethod
	def _read_bam(path):
		save = pysam.set_verbosity(0)
		f1 = pysam.AlignmentFile(path, 'rb', check_sq=False)
		pysam.set_verbosity(save)
		for line in f1:
			yield line.query_name, line.query_sequence, line.query_qualities, None
		f1.close()

	@staticmethod
	def _read_seq_data(path):
		if path.endswith('.fastq.gz') or path.endswith('.fq.gz') or path.endswith('.fastq') or path.endswith('.fq'):
			for rid, seq, qual, info in QC_IO._read_fq(path):
				yield rid, seq, qual, info
		elif path.endswith('.bam'):
			for rid, seq, qual, info in QC_IO._read_bam(path):
				yield rid, seq, qual, info

	@staticmethod
	def mkdir(out_dir):
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

	@staticmethod
	def rm_info_files(out_dir):
		if not os.path.exists(out_dir):
			return
		info_files = os.listdir(out_dir)
		for x in info_files:
			if x.endswith('.info.txt'):
				os.unlink(os.path.join(out_dir, x))


class NpVec(object):

	def __init__(self, init_size, init_ele, dtype):
		if init_size <= 0:
			self._nparr = np.array([], dtype=dtype)
		else:
			self._nparr = np.array([init_ele], dtype=dtype)
			self._nparr = np.resize(self._nparr, init_size)
		self._real_size = 0
		self._init_ele = init_ele
		self._dtype = dtype
	
	def append(self, ele):
		if self._real_size >= len(self._nparr):
			new_size = int(self._real_size * 1.2 + 0.5) + 1 if self._real_size > 0 else 4
			self._nparr = np.resize(self._nparr, new_size)
			gc.collect()
			for i in range(self._real_size, new_size):
				self._nparr[i] = self._init_ele
		self._nparr[self._real_size] = ele
		self._real_size += 1

	def shrink(self):
		self._nparr = np.resize(self._nparr, self._real_size)
		gc.collect()

	def clear(self):
		self._nparr = np.array([], dtype=self._dtype)
		gc.collect()
		self._real_size = 0


# su_infos: read_id, qscore, pass_filter
def read_write_summary(su_paths, su_infos, LOG):
	COLNAMES = ['read_id', 'mean_qscore_template', 'passes_filtering']
	out_header = None
	rid_idx, qs_idx, pass_idx = None, None, None
	for path in su_paths:
		with read_file(path) as fin:
			for line in fin:
				if line.startswith('filename') or line.startswith('read'):
					header = line.strip().split()
					if out_header is not None and len(header) != len(out_header):
						raise Exception('the header is not the same as before in %s' % path)
					for colname in COLNAMES:
						if colname not in header:
							raise Exception('the column %s is not found in %s' % (colname, ','.join(header)))
					if out_header is None:
						rid_idx, qs_idx, pass_idx = header.index(COLNAMES[0]), header.index(COLNAMES[1]), header.index(COLNAMES[2])
						out_header = header
					else:
						new_rid_idx, new_qs_idx, new_pass_idx = header.index(COLNAMES[0]), header.index(COLNAMES[1]), header.index(COLNAMES[2])
						if new_rid_idx != rid_idx or new_qs_idx != qs_idx or new_pass_idx != pass_idx:
							raise Exception('the header is not the same as before in %s' % path)
					continue
				field = line.strip().split()
				if len(field) != len(out_header):
					LOG.warning('the colnum is not equal with header in %s' % path)
					continue
				read_id = field[rid_idx]
				read_id = transform_read_id(read_id)
				qscore = float(field[qs_idx])
				pass_filter = True if field[pass_idx] == "TRUE" else False
				su_infos.append((read_id, qscore, pass_filter))
	LOG.info('Done reading summary file')


def transform_read_id(read_id):
	return read_id.replace('-', '')


def get_su_i(np_su_infos, read_id):
	if read_id < np_su_infos[0][0]:
		return None
	if read_id > np_su_infos[-1][0]:
		return None
	s, e = 0, len(np_su_infos) - 1
	while s <= e:
		m = (s + e) // 2
		if read_id == np_su_infos[m][0]:
			return m
		elif read_id < np_su_infos[m][0]:
			e = m - 1
		else:
			s = m + 1


def get_prefix(fq_file):
	prefix = os.path.basename(fq_file)
	prefix = re.sub('.fastq.gz$', '', prefix)
	prefix = re.sub('.fastq$', '', prefix)
	prefix = re.sub('.fq.gz$', '', prefix)
	prefix = re.sub('.fq$', '', prefix)
	prefix = re.sub('.bam$', '', prefix)
	return prefix


def normal_run(args, LOG):
	QC_IO.rm_info_files(args.out_dir)
	QC_IO.mkdir(args.out_dir)

	su_paths = QC_IO.read_lst_file(args.summary_lst)
	dt = np.dtype([('read_id', 'S32'), ('qscore', np.float32), ('pass_filter', np.bool_)])
	init_ele = ('', 0., True)
	su_infos = NpVec(10000, init_ele, dt)  # read_id, qscore, pass_filter
	read_write_summary(su_paths, su_infos, LOG)
	su_infos.shrink()
	np_su_infos = su_infos._nparr
	np_su_infos.sort(order='read_id')

	fq_paths = QC_IO.read_lst_file(args.fq_lst)
	for i, fq_file in enumerate(fq_paths):
		prefix = '%s.%s' % (get_prefix(fq_file), i)
		out_info_prefix = os.path.join(args.out_dir, '%s.info.txt' % prefix)
		with write_file(out_info_prefix) as fout:
			for rid, seq, qual, info in QC_IO._read_seq_data(fq_file):
				gc_count = seq.upper().count('G') + seq.upper().count('C')
				rid = transform_read_id(rid)
				j = get_su_i(np_su_infos, to_byte(rid))
				if j is None:
					LOG.warning('no summary info for read id %s' % rid)
					continue
				# read_id, gc_count, length, qscore, pass_filter
				fout.write('{"read_id": "%s", "gc_count": %s, "length": %s, "qscore": %.6f, "pass_filter": "%s"}\n' %
					(rid, gc_count, len(seq), float(np_su_infos[j][1]), "TRUE" if np_su_infos[j][2] else "FALSE"))
		LOG.info('Done dealing with fastq/bam file %s' % fq_file)
	np_su_infos = None
	su_infos.clear()
	gc.collect()


class StatAndFig():

	@staticmethod
	def plot_gc(plt, _gc_percent_list, out_img_file, bins=100):
		fig, ax = plt.subplots(figsize=(10, 6), )
		ax.hist(_gc_percent_list, bins=bins, range=(0, 100), density=False, alpha=0.75, edgecolor='white',
				linewidth=0.5)
		ax.set_xlabel('GC_content')
		ax.set_ylabel('Read Count')
		plt.xlim(0, 100)
		plt.savefig(out_img_file, dpi=700)
		plt.close()

	@staticmethod
	def plot_reads_dis(plt, _length_list, out_img_file, xmin=0, xmax=300 * 1000, bins=4000):
		fig, ax = plt.subplots(figsize=(10, 6), )
		ax.hist(_length_list, bins, density=False, alpha=0.75, )
		ax.set_xlabel('Read Length')
		ax.set_ylabel('Read Count')
		plt.xlim(xmin, xmax)
		plt.savefig(out_img_file, dpi=700)
		plt.close()

	@staticmethod
	def plot_qscore_length_dis(plt, transforms, _length_list, _qscore_list, out_img_file):
		# generate counts
		len_part_size, q_part_size = 5000, 10
		if len(_length_list) > 0:
			minl, maxl = 0, int((max(_length_list) - 1) / len_part_size) + 1
		else:
			minl, maxl = 0, 1
		if len(_qscore_list) > 0:
			minq, maxq = 0, int((max(_qscore_list) - 0.000001) * q_part_size) + 1
		else:
			minq, maxq = 0, 1
		# print(minl, maxl, minq, maxq)
		counts = defaultdict(dict)
		num = min(len(_length_list), len(_qscore_list))
		lcounts = defaultdict(int)
		qcounts = defaultdict(int)
		for i in range(num):
			lval = int(_length_list[i] / len_part_size)
			qval = int(_qscore_list[i] * q_part_size)
			if lval > maxl:
				lval = maxl
				#LOG.info('length, %s %s %s' % (lval, maxl, _length_list[i]))
			if qval > maxq:
				qval = maxq
				#LOG.info('qscore, %s %s %s' % (qval, maxq, _qscore_list[i]))
			if qval not in counts[lval]:
				counts[lval][qval] = 1
			else:
				counts[lval][qval] += 1
			lcounts[lval] += 1
			qcounts[qval] += 1

		# create data for figure
		X = [x * len_part_size for x in range(minl, maxl + 1)]
		Y = [y / q_part_size for y in range(minq, maxq + 1)]
		Z = [[0 for y in range(minl, maxl + 1)] for x in range(minq, maxq + 1)]
		for lval, vals in counts.items():
			for qval, count in vals.items():
				Z[qval][lval] = count
		lz = [0 for y in range(minl, maxl + 1)]
		for lval, count in lcounts.items():
			lz[lval] = count
		qz = [0 for x in range(minq, maxq + 1)]
		for qval, count in qcounts.items():
			qz[qval] = -count
		# plt.style.use('ggplot')

		# mid plot
		fig = plt.figure(figsize=(5,5))
		plt.rcParams['axes.facecolor'] = '#EAEAF2'
		plt.rcParams['axes.labelsize'] = 7
		# plt.rcParams['font.size'] = 9
		mid = fig.add_axes([0.1, 0.1, 0.7, 0.7], xlim=(minl, maxl * len_part_size), ylim=(minq, maxq / q_part_size))
		mid.grid(True, color='white', linestyle='-', linewidth=0.5, alpha=0.8, zorder=0)
		mid.patch.set_facecolor('#EAEAF2')
		# cs = mid.contourf(X, Y, Z, levels=8, cmap='Greens', extend='both', alpha=0.8)
		curr_cmap = copy.copy(plt.cm.get_cmap("Greens"))
		cs = mid.contourf(X, Y, Z, levels=8, cmap=curr_cmap, extend='both', alpha=0.8)
		cs.cmap.set_under('#EAEAF2')		
		cs.changed()
		mid.set_xlabel('Read lengths')
		mid.set_ylabel('Average read quality scores')
		mid.spines['top'].set_visible(False)
		mid.spines['right'].set_visible(False)
		mid.spines['bottom'].set_visible(False)
		mid.spines['left'].set_visible(False)
		mid.tick_params(color='white', labelsize=7)
		# C = mid.contour(X, Y, Z, colors='blue')
		# mid.clabel(C, inline=1, fontsize=5)

		# top plot
		max_lz = max(lz)
		if max_lz == 0: max_lz = 1
		top = fig.add_axes([0.1, 0.805, 0.7, 0.13], xlim=(minl, maxl * len_part_size), ylim=(0, max_lz*1.2))
		# top.set_title('Average read quality vs Read lengths plot')
		top.text(0.6, 1.1, 'Average read quality vs Read lengths plot', transform=top.transAxes, ha='center')
		top.grid(True, color='white', linestyle='-', linewidth=0.5, alpha=0.8, zorder=0)
		top.patch.set_facecolor('#EAEAF2')
		top.fill_between(X, lz, color='#C2ECD9', alpha=0.8)
		top.plot(X, lz, '-c', linewidth=1)
		top.set_xticklabels([])
		top.set_yticklabels([])
		top.spines['top'].set_visible(False)
		top.spines['right'].set_visible(False)
		top.spines['bottom'].set_visible(False)
		top.spines['left'].set_visible(False)
		top.tick_params(color='white')

		# right plot
		min_qz = min(qz)
		if min_qz == 0: min_qz = 1
		right = fig.add_axes([0.805, 0.1, 0.13, 0.7], ylim=(minq, maxq / q_part_size), xlim=(0, -min_qz*1.2))
		right.grid(True, color='white', linestyle='-', linewidth=0.5, alpha=0.8, zorder=0)
		right.patch.set_facecolor('#EAEAF2')
		base = right.transData
		rot = transforms.Affine2D().rotate_deg(90)
		right.fill_between(Y, qz, color='#C2ECD9', transform=rot + base, alpha=0.8)
		right.plot(Y, qz, '-c', linewidth=1, transform=rot + base)
		right.set_yticklabels([])
		right.set_xticklabels([])
		right.spines['top'].set_visible(False)
		right.spines['right'].set_visible(False)
		right.spines['bottom'].set_visible(False)
		right.spines['left'].set_visible(False)
		right.tick_params(color='white')

		plt.savefig(out_img_file, dpi=700)
		plt.close()

	@staticmethod
	def comp_nx0_on_sorted_lens(sorted_lengths, base_num, perc=0.5):
		if len(sorted_lengths) == 0:
			return 0, 0
		half_length = perc * base_num
		nx0 = 0
		sum_len = 0
		count = 0
		for length in sorted_lengths:
			sum_len += length
			count += 1
			if sum_len >= half_length:
				nx0 = length
				break
		return nx0, count

	# input is reversely sorted length list
	@staticmethod
	def nx0_v2(_length_list, base_num):
		_nx0 = dict()
		for level_num in range(1, 10):
			level = 'n%d0' % level_num
			nx0_val = int(StatAndFig.comp_nx0_on_sorted_lens(_length_list, base_num, level_num / 10)[0])
			if nx0_val == 0: nx0_val = 1
			_nx0[level] = nx0_val
		return _nx0

	@staticmethod
	def nx0_v1(_length_list, _base_num):
		_nx0 = dict()
		for level_num in range(1, 10):
			level = 'n%d0' % level_num
			_nx0[level] = 0
		tmp_base_num = 0
		for length in _length_list:
			tmp_base_num += length
			level_num = int(10 * tmp_base_num / _base_num)
			if level_num == 0 or level_num == 10: continue
			level = 'n%d0' % level_num
			if _nx0.get(level) == 0: _nx0[level] = int(length)
		return _nx0
	
	@staticmethod
	def nx0(_length_list, base_num):
		if len(_length_list) > 10000:
			return StatAndFig.nx0_v1(_length_list, base_num)
		else:
			return StatAndFig.nx0_v2(_length_list, base_num)

	@staticmethod
	def ultralong(_length_list):
		_ultralong = dict()
		# >50kb，>100kb，>150kb，>200kb，>300kb，>400kb，>500kb的产量及reads数
		data = dict()
		for length in _length_list:
			key = int((length - 1) / 10000)
			data.setdefault(key, [0, 0])
			data[key][0] += int(length)
			data[key][1] += 1
		for threshold in [5, 10,12, 15, 20, 30, 40, 50]:
			level = '>%d0kb' % threshold
			_ultralong.setdefault(level, {'base_num':0, 'read_num':0})
			for key, value in data.items():
				if key >= threshold:
					_ultralong[level]['base_num'] += value[0]
					_ultralong[level]['read_num'] += value[1]
		return _ultralong

	# _length_list is reversely sorted
	@staticmethod
	def n50_gt(_length_list):
		_n50_gt = dict()
		a = [5, 10, 12, 15, 20]
		for a_i in range(len(a)):
			lt_base_num = gt_base_num = 0
			for i in _length_list:
				if (i >= a[a_i] * 10000):
					gt_base_num += i
				else:
					lt_base_num += i
					if (lt_base_num >= gt_base_num):
						break
			_n50_gt["n50>%d0kb" % a[a_i]] = int(lt_base_num + gt_base_num)
		return _n50_gt

	# warning: infos_wo_rid[1] are modified which are sorted reversely
	@staticmethod
	def write_json(infos_wo_rid, out_qc_stat_json_file):
		# import pdb; pdb.set_trace()
		pass_lens = np.compress(infos_wo_rid[3], infos_wo_rid[1])
		pass_lens.sort()
		pass_lens = pass_lens[::-1]
		pass_base_num = int(np.sum(pass_lens))
		pass_nx0 = StatAndFig.nx0(pass_lens, pass_base_num)
		pass_ultralong = StatAndFig.ultralong(pass_lens)
		pass_n50_gt = StatAndFig.n50_gt(pass_lens)
		pass_read_num = len(pass_lens)
		pass_max_len = int(np.max(pass_lens)) if len(pass_lens) > 0 else 1
		pass_mean_len = float(np.mean(pass_lens)) if len(pass_lens) > 0 else 1.
		pass_medium_len = float(np.median(pass_lens)) if len(pass_lens) > 0 else 1.
		pass_lens = None; gc.collect()

		pass_qscore = np.compress(infos_wo_rid[3], infos_wo_rid[2])
		pass_mean_qscore = float(np.mean(pass_qscore)) if len(pass_qscore) > 0 else 1.
		pass_qscore = None; gc.collect()

		pass_gc = np.compress(infos_wo_rid[3], infos_wo_rid[0])
		pass_gc_rate = float(np.sum(pass_gc) / pass_base_num) if pass_base_num > 0 else 1.
		pass_gc = None; gc.collect()

		all_lens = infos_wo_rid[1]
		all_lens.sort()
		all_lens = all_lens[::-1]
		all_base_num = int(np.sum(all_lens))
		total_read_num = len(all_lens)
		all_nx0 = StatAndFig.nx0(all_lens, all_base_num)
		all_lens = None; gc.collect()

		fo = open(out_qc_stat_json_file, 'w')
		fo.write(json.dumps({
			'total_base_num': all_base_num if all_base_num > 0 else 1,
			'total_read_num': total_read_num if total_read_num > 0 else 1,
			'pass_base_num': pass_base_num if pass_base_num > 0 else 1,
			'pass_read_num': pass_read_num if pass_read_num > 0 else 1,
			'pass_max_len': pass_max_len,
			'pass_mean_len': pass_mean_len,
			'pass_medium_len': pass_medium_len,
			'pass_mean_qscore': pass_mean_qscore,
			'pass_n50_len': pass_nx0['n50'],
			'pass_gc_rate': pass_gc_rate,
			'all_nx0': all_nx0,
			'pass_nx0': pass_nx0,
			'pass_ultralong': pass_ultralong,
			'pass_n50>x_base_num': pass_n50_gt
		}, indent=4))
		fo.close()

	@staticmethod
	def read_info(info_file_lst):
		info_files = set([os.path.abspath(x) for x in QC_IO.read_lst_file(info_file_lst)])
		infos_wo_rid = [array.array('l'), array.array('l'), array.array('d'), array.array('b')]  # gc_count, length, qscore, pass_filter
		for info_file in info_files:
			with read_file(info_file) as fin:
				for line in fin:
					line = line.strip()
					if not line: continue
					info = json.loads(line)
					infos_wo_rid[0].append(info['gc_count'])
					infos_wo_rid[1].append(info['length'])
					infos_wo_rid[2].append(info['qscore'])
					infos_wo_rid[3].append(1 if info['pass_filter'] == 'TRUE' else 0)
		ret_arr = [np.array(infos_wo_rid[0], dtype=np.int32), np.array(infos_wo_rid[1], dtype=np.int32),
			np.array(infos_wo_rid[2], dtype=np.float64), np.array(infos_wo_rid[3], dtype=np.bool_)]
		infos_wo_rid.clear()
		gc.collect()
		return ret_arr

	@staticmethod
	def create_stats_and_figs(info_file_lst, out_stat_dir):
		QC_IO.mkdir(out_stat_dir)
		infos_wo_rid = StatAndFig.read_info(info_file_lst)
		StatAndFig.write_json(infos_wo_rid, os.path.join(out_stat_dir, 'qc_stat.json'))

		import matplotlib
		matplotlib.use('Agg')
		from matplotlib import pyplot as plt, transforms
		out_gc_content_img_file = os.path.join(out_stat_dir, 'GC_content.png')
		out_reads_dis_img_file = os.path.join(out_stat_dir, 'reads_distribution.png')
		out_qscore_length_dis_img_file = os.path.join(out_stat_dir, 'qscore_vs_length.png')

		pass_lens = np.compress(infos_wo_rid[3], infos_wo_rid[1])
		pass_gc = np.compress(infos_wo_rid[3], infos_wo_rid[0])
		gc_percent = pass_gc / pass_lens
		pass_gc = None; gc.collect()
		StatAndFig.plot_gc(plt, gc_percent, out_gc_content_img_file)
		gc_percent = None; gc.collect()
		StatAndFig.plot_reads_dis(plt, pass_lens, out_reads_dis_img_file)
		pass_qscore = np.compress(infos_wo_rid[3], infos_wo_rid[2])
		StatAndFig.plot_qscore_length_dis(plt, transforms, pass_lens, pass_qscore, out_qscore_length_dis_img_file)
		pass_qscore = pass_lens = None; gc.collect()


def main(args):
	logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(asctime)s %(message)s')
	LOG = logging.getLogger(__name__)
	LOG.info(args)
	try:
		normal_run(args, LOG)
	except:
		LOG.error(traceback.format_exc())
		raise


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='GrandOmics Copyright(c) 2022. Author: Zhuo Wang;'
	)
	parser.add_argument('fq_lst', help='input fastq/bam path list')
	parser.add_argument('summary_lst', help='input summary path list')
	parser.add_argument('out_dir', help='output directory')
	parser.set_defaults(function=main)
	args = parser.parse_args()
	args.function(args)
