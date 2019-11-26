# Question 1-5： BAM FLAG
## 1. Explain BAM FLAG value： 143
+ 1: template having multiple segments in sequencing
+ 2: each segment properly aligned according to the aligner
+ 4: segment unmapped
+ 8: next segment in the template unmapped
+ 128: the last segment in the template

## 2. Explain BAM FLAG value： 99
+ 1: template having multiple segments in sequencing
+ 2: each segment properly aligned according to the aligner
+ 32: SEQ of the next segment in the template being reverse complemented
+ 64: the first segment in the template

## 3. Explain BAM FLAG value：516
+ 4: segment unmapped
+ 512: not passing filters, such as platform/vendor quality controls

## 4. Explain BAM FLAG value： 2064
+ 16: SEQ being reverse complemented
+ 2048: supplementary alignment

## 5. Explain BAM FLAG value： 147
+ 1: template having multiple segments in sequencing
+ 2: each segment properly aligned according to the aligner
+ 16: SEQ being reverse complemented
+ 128: the last segment in the template
