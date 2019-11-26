#encoding='utf-8'
import os,sys
import argparse

bitInfoArr = [
    'template having multiple segments in sequencing--()',
    'each segment properly aligned according to the aligner',
    'segment unmapped',
    'next segment in the template unmapped',
    'SEQ being reverse complemented',
    'SEQ of the next segment in the template being reverse complemented',
    'the first segment in the template',
    'the last segment in the template',
    'secondary alignment',
    'not passing filters, such as platform/vendor quality controls',
    'PCR or optical duplicate',
    'supplementary alignment'
]

def getCigarInfo(flag):
    for i in range(0, len(bitInfoArr)):
        if flag & (1<<i) > 0:
            print(str(1<<i)+': '+bitInfoArr[i])

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
