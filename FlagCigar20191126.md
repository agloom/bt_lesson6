# Question 1-5： BAM FLAG
## 1. Explain BAM FLAG value： 143
+ 1: template having multiple segments in sequencing---测序片段不止一条（PE测序）
+ 2: each segment properly aligned according to the aligner---和参考序列完全匹配
+ 4: segment unmapped---序列没有比对上参考序列上
+ 8: next segment in the template unmapped---另一端片段没有比对到参考序列上
+ 128: the last segment in the template---当前序列测序序列的read2

## 2. Explain BAM FLAG value： 99
+ 1: template having multiple segments in sequencing---测序片段不止一条（PE测序）
+ 2: each segment properly aligned according to the aligner---和参考序列完全匹配
+ 32: SEQ of the next segment in the template being reverse complemented---另一端序列比对到参考序列负链上
+ 64: the first segment in the template---当前序列为测序序列的read1

## 3. Explain BAM FLAG value：516
+ 4: segment unmapped---序列没有比对上参考序列上
+ 512: not passing filters, such as platform/vendor quality controls---过滤标志，表示这个序列质量较差

## 4. Explain BAM FLAG value： 2064
+ 16: SEQ being reverse complemented---序列比对到参考序列负链
+ 2048: supplementary alignment---序列补充比对

## 5. Explain BAM FLAG value： 147
+ 1: template having multiple segments in sequencing---测序片段不止一条（PE测序）
+ 2: each segment properly aligned according to the aligner---和参考序列完全匹配
+ 16: SEQ being reverse complemented---序列比对到参考序列负链
+ 128: the last segment in the template---当前序列测序序列的read2

# Question 6-10： BAM CIGAR
## 6. Explain BAM CIGAR：14M2D31M
+ 按序列顺序，14个碱基匹配2个碱基缺失31个碱基匹配

## 7. Explain BAM CIGAR：3S6M1D5M
+ 3个碱基属于soft clipping，接着6个碱基匹配1个碱基缺失5个碱基匹配

## 8. Explain BAM CIGAR: 6M14N5M
+ 6个碱基匹配14个碱基跳过（skipped, 对于mRNA，为内含子）5个碱基匹配

## 9. Explain BAM CIGAR: 7M5D8M2I14M  (小写：7m5d8m2i14m）
+ 7个碱基匹配5个碱基缺失8个碱基匹配2个碱基插入14个碱基匹配

## 10. how long is the read with alignment CIGAR of 7M5D8M2I14M?
31bp (7+8+2+14)
