#encoding: utf-8
#author: zhenweibo
#email: agloom@163.com
import os,sys
import argparse

bitInfoArr = [
    'template having multiple segments in sequencing---测序片段不止一条（PE测序）',
    'each segment properly aligned according to the aligner---和参考序列完全匹配',
    'segment unmapped---序列没有比对上参考序列上',
    'next segment in the template unmapped---另一端片段没有比对到参考序列上',
    'SEQ being reverse complemented---序列比对到参考序列负链',
    'SEQ of the next segment in the template being reverse complemented---另一端序列比对到参考序列负链上',
    'the first segment in the template---当前序列为测序序列的read1',
    'the last segment in the template---当前序列测序序列的read2',
    'secondary alignment---次要比对位置，在序列比对到参考序列多个位置上会出现。区别于主要比对位置',
    'not passing filters, such as platform/vendor quality controls---过滤标志，表示这个序列质量较差',
    'PCR or optical duplicate---重复的序列，一般用于picard等软件用于标记由于PCR等导致的重复序列',
    'supplementary alignment---序列补充比对'
]

def getCigarInfo(flag):
    for i in range(0, len(bitInfoArr)):
        if flag & (1<<i) > 0:
            print('+ '+str(1<<i)+': '+bitInfoArr[i])

def main():
    parser = argparse.ArgumentParser(description='get bam cigar information.')
    #parser.add_argument('-c', '--cigar', dest='cigar', required=True, help='cigar')
    parser.add_argument('-f', '--flag', dest='flag', required=True, help='flag number')
    argv = parser.parse_args()
    #cigar = argv.cigar.strip()
    flag = int(argv.flag.strip())
    getCigarInfo(flag)

if __name__=='__main__':
    main()
