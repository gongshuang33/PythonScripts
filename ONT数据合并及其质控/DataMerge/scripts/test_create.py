import sys
sys.path.append('/home/dm/nextdataflow3/ONT/qc/')
from ont_qc_info import StatAndFig
#import logging

#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(asctime)s %(message)s')
#LOG = logging.getLogger(__name__)

info_file_lst = 'info.lst'
#info_file_lst = 'null.lst'
out_stat_dir = './'
StatAndFig.create_stats_and_figs(info_file_lst, out_stat_dir)

