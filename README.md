# bt_lesson6
job for bioinformatic lesson6

# SyntenyPlot

	usage:  nqc -h show the usage detail.

## Descriptions:

	nqc is a multi-threads program, written by c program language, for ngs fastq data quality control. 

## Requirements: 

	gcc supported C99

#HOW-TO

### 1. **compile**
gcc -o nqc main.c -lz

### 2. **usage**
You can execute nqc -h to get usage information.
    Usage:
     -f    the path of fq files or fq.gz files, seperated by comma.
     
     -p     threads number.
     
     -o     output directory.
    
     -s     sample name.
     
     -N     N rate. default: 0.1"
     
     -l     low qual.
     
     -L     the minimal length of read. The reads(perhaps have been truncated) will be discarded if its length is less than this value.
     
     -g     the sequence data is compressed by gzip.
     
     -h    help information.


### TO-DO list;

+ [ ] 1. trim adaptors


## Citation:

No

## Author:

---------------------------------------------------------------------

>	**zhenweibo**
>	E-mail: <agloom@163.com>

---------------------------------------------------------------------

## Copyright

Copyright (c) 2019-2021 zhenweibo

Under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version, Permission is hereby granted, 
free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, 
copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
