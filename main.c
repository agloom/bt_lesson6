/**
author: zhenweibo
email: agloom@163.com
date: 2016-08-21
description: quality control program for NGS data.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>             //std=c99编译时需要加这一个头文件
#include <pthread.h>
#include <semaphore.h>
#include "zlib.h"
#include "qc.h"

#define PATH_MAX_LEN  2048     //路径最长字节数
#define NAME_LEN    256             //输出名称前缀长度
const unsigned int INIT_READ_LEN = 1024 ;      //在不知道reads长度的情况下，用于初始长度值
const unsigned int EACH_READ_NUM = 20000 ;  //每个线程每次处理的reads数目。

unsigned int FQ_NUM = 2 ;            //初始化为fq文件最大数目，后面在主线程中会根据实际的文件数目更改此变量(不超过初始值)
unsigned int READ_LEN = 0;      //读取fq文件后，获取的reads长度(包含最后面一个换行符)保存到此变量中
unsigned int LINE_LEN = 0;      //read行的长度，包含后面换行符和\0
unsigned int EACH_LINE_NUM = 0;         //每个read占4行，在main函数中再赋值
int PHRED_BASE = 0;                             //phred的基数
unsigned int READ_MIN_BASES = 0;            //如果截短reads，剩下的长度如果小于这个值，则会丢弃整个reads
float   N_Rate = 0.1;                           //每个read里允许最多的N比例。
unsigned int LOW_QUAL = 5 ;             //低质量值
char out_dir[PATH_MAX_LEN] ;            //输出目录
char sample[NAME_LEN];                //样本名称
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
unsigned int finished_thread_num = 0;

/**
程序说明
*/
void usage(){
    puts("***************");
    puts("Usage: ");
    puts("\t -f    the path of fq files or fq.gz files, seperated by comma.");
    puts("\t -p     threads number.");
    puts("\t -o     output directory.");
    puts("\t -s     sample name.");
    puts("\t -N     N rate. default: 0.1");
    puts("\t -l     low qual.");
    puts("\t -L     the minimal length of read. The reads(perhaps have been truncated) will be discarded if its length is less than this value.");
    puts("\t -g     the sequence data is compressed by gzip.");
    puts("\t -h    help information.");
    puts("***************");
}

/**
根据fq内容判断phred33还是phred64
*/
int phredVal(char *path, unsigned int rlen){
    int len = rlen > 60 ? rlen : 60;
    len = len +2 ;    // '\n\0'
    char  buf[len];
    char *rb;
    gzFile fqFile = gzopen(path, "r");
    while(!gzeof(fqFile)){
        memset(buf, 0, len);
        for(int i=0; i<4; i++){
            rb = gzgets(fqFile, buf, len);
            if(rb == NULL){
                fprintf(stderr, "Failed to get line from %s.\n", path);
                gzclose(fqFile);
                return -1;
            }
            //如果读取了第4行（质量值）
            if(i==3){
                for(int n =0; n< rlen; n++){
                    if(buf[n]< 64 && buf[n]>32){
                        return 33;
                    }else if(buf[n] > 74 && buf[n]<105){
                        return 64;
                    }else if(buf[n]<33 || buf[n] >104){
                        return -1;
                    }
                }
            }
        }
    }
    return -1;
}

/**
判断数据中所有的文件流是否已达到末尾。只要有一个未到达，则返回false
*/
_Bool isAllFileEnd(gzFile  a[], unsigned int n){
    for(int i=0; i<n; i++){
        if(!gzeof(a[i])){
            return  false;
        }
    }
    return true;
}
/***** 主函数 *****/
int main(int argc, char** argv)
{
    EACH_LINE_NUM = 4*EACH_READ_NUM;         //每个read占4行
    char fqs[PATH_MAX_LEN];              //输入的fq文件路径，多个fq路径之间以英文逗号分隔
    memset(fqs, 0, PATH_MAX_LEN);
    memset(out_dir, 0, PATH_MAX_LEN);

    int thread_num = 5;                     //thread number
    _Bool bGz = false;                       //输入的fq文件是否为gz压缩格式

    //get parameters
    char opt;
    while((opt = getopt(argc, argv, "f:L:l:s:N:o:gp:h")) != -1){
        switch(opt){
            case 'f':
                if(strlen(optarg) > PATH_MAX_LEN-1){
                    puts("The length of fqs file paths is too long.\n");
                    exit(2);
                }

                strncpy(fqs, optarg, PATH_MAX_LEN);
                break;
            case 'L':
                READ_MIN_BASES = atoi(optarg);
                if(READ_MIN_BASES <0){
                    fprintf(stderr, "-L can not be less than 0.\n");
                    exit(2);
                }
                break;
            case 's':
                memcpy(sample, optarg, strlen(optarg)+1<NAME_LEN? strlen(optarg)+1:NAME_LEN);
                break;
            case 'l':
                LOW_QUAL = atoi(optarg);
                break;
            case 'o':
                if(strlen(optarg)+1 > PATH_MAX_LEN){
                    fprintf(stderr, "The output path is too long.\n");
                    exit(2);
                }
                strncpy(out_dir,optarg, strlen(optarg)+1);
                break;
            case 'N':
                N_Rate = atof(optarg);
                if(N_Rate<0 || N_Rate>1){
                    fprintf(stderr, "Error: -N must be in 0~1.\n");
                    exit(2);
                }
                break;
            case 'g':
                bGz = true;
                break;
            case 'p':
                thread_num = atoi(optarg);
                if(thread_num <=0){
                    fprintf(stderr, "Warning: the value of -t is invalid. use default value.\n");
                }
                break;
            case 'h':
                usage();
                exit(-1);
            default:
                usage();
                exit(-1);
        }
    }

    //*************************
    char fqA[FQ_NUM][PATH_MAX_LEN];              //存放fq路径的数组
    int fqNum = 0;                                              //实际fq文件数目
    if(access(out_dir, F_OK|W_OK)!=0){
        fprintf(stderr, "Error: output directory %s should exists and can be written.\n", out_dir);
        exit(2);
    }
    //*************************

    for(int i =0; i< FQ_NUM; i++){
        memset(fqA[i], 0, PATH_MAX_LEN);
    }

    char *fqPath;
    char *fqp = fqs;
    while((fqPath = strtok(fqp, ",")) != NULL && fqNum < FQ_NUM){
        strncpy(fqA[fqNum], fqPath, PATH_MAX_LEN);
        ++fqNum ;
        fqp = NULL;
    }
    //更改FQ_NUM为实际的fq文件数量
    FQ_NUM = fqNum;

    /*******************
    计算read长度
    *******************/

    gzFile GFILE = gzopen(fqA[0], "r");
    if(GFILE == NULL){
        fprintf(stderr , "Error! Failed to open %s\n", fqA[0]);
        exit(3);
    }
    unsigned int lineNo = 1;
    while(!gzeof(GFILE)){
        char *rb;
        char line[INIT_READ_LEN];
        memset(line, 0, INIT_READ_LEN);
        rb = gzgets(GFILE, line, INIT_READ_LEN);
        if(rb == NULL){
            fprintf(stderr, "Error: can not get the first line content of %s\n.", fqA[0]);
            exit(3);
        }
        if(lineNo == 2){        //从第二行开始计算read长度
            READ_LEN += strlen(line);
        }
//        printf("tmp_len=%d\n", READ_LEN);
        //如果已经读取完整的一行
        if(line[INIT_READ_LEN-2] == '\0' || line[INIT_READ_LEN-2] == '\n'){
            lineNo++;
            if(lineNo > 2){     //第二行结束s
                break;
            }
        }
    }

    if(READ_LEN == 0){
        fprintf(stderr, "Error ! Failed to get read length.\n");
        exit(3);
    }
    LINE_LEN = READ_LEN +1;
    gzclose(GFILE);
    printf("read length: %d\n", READ_LEN-1);

    /******判断为phred33还是phred64*****/
    PHRED_BASE = phredVal(fqA[0], READ_LEN-1);
    if(PHRED_BASE !=33 && PHRED_BASE !=64){
        fprintf(stderr, "Unknown phred: %d\n", PHRED_BASE);
        exit(7);
    }

    /*******************
    创建线程
    *******************/

    /****** 存放线程输入输出的数据缓存 *****/
    char *pInBufArr[thread_num][fqNum];     //centOS6.7 64位系统中int 4字节，指针类型8字节
    char *pOutBufArr[thread_num][fqNum];
    memset(pInBufArr, 0 , sizeof(pInBufArr));
    memset(pOutBufArr, 0 , sizeof(pInBufArr));
    //给每个指针分配内存
    for(int i=0; i< thread_num; i++){
        for(int j =0; j< fqNum; j++){
            pInBufArr[i][j] = malloc(EACH_LINE_NUM*(LINE_LEN)*sizeof(char));
            memset(pInBufArr[i][j], 0, EACH_LINE_NUM*(LINE_LEN)*sizeof(char));
            pOutBufArr[i][j] = malloc(EACH_LINE_NUM*(LINE_LEN)*sizeof(char));
            memset(pOutBufArr[i][j], 0, EACH_LINE_NUM*(LINE_LEN)*sizeof(char));

            if(pInBufArr[i][j] == NULL){
                fprintf(stderr, "Error! Failed to allocate memory(input) to thread%d\n", i);
                exit(4);
            }
            if(pOutBufArr[i][j] == NULL){
                fprintf(stderr, "Error! Failed to allocate memory(output) to thread%d\n", i);
                exit(4);
            }
        }
    }

    pthread_t  th[thread_num];
    unsigned int t_flag[thread_num];                            //每个线程的运行标志
    sem_t inSemEmpty[thread_num];            //信号量，用于主线程与子线程(处理线程)之间
    sem_t inSemFull[thread_num];
    sem_t outSemEmpty[thread_num];           //信号量，用于保存线程与子线程(处理线程)之间
    sem_t outSemFull[thread_num];

    TParam param[thread_num];       //存放各线程参数
    TResult result[thread_num];         //存放线程结果
    unsigned long long posBase[thread_num][(READ_LEN-1)*FQ_NUM][5];  //[线程][位置][碱基类别atgcn]
    unsigned long long posPhred[thread_num][(READ_LEN-1)*FQ_NUM];
    for (int i = 0; i< thread_num; i++){
        //初始化posPhred和posBase
        for(int j=0; j<(READ_LEN-1)*FQ_NUM; j++){
            posPhred[i][j] = 0;

            for(int x=0; x<5; x++){
                posBase[i][j][x] = 0;
            }
        }
        result[i].posPhred = posPhred[i];
        result[i].posBase = posBase[i];
        //初始化运行状态标志
        t_flag[i] = 1;

        //-------------打印信号量地址--------------
        printf("########\nThread%d Sem(from main):\n inEmpty: %lu\tinFull: %lu\toutEmpty: %lu\toutFull: %lu\n#########\n", i, &inSemEmpty[i], &inSemFull[i], &outSemEmpty[i], &outSemFull[i]);
        //---------------------------

        //初始化信号量
        sem_init(&inSemEmpty[i], 0, 1);
        sem_init(&inSemFull[i], 0, 0);
        sem_init(&outSemEmpty[i], 0, 1);
        sem_init(&outSemFull[i], 0, 0);

        param[i].cId = i;
        param[i].pIn = pInBufArr[i];
        param[i].pOut = pOutBufArr[i];
        param[i].inSemEmpty = inSemEmpty;
        param[i].inSemFull = inSemFull;
        param[i].outSemEmpty = outSemEmpty;
        param[i].outSemFull = outSemFull;
        param[i].runFlag = &t_flag[i];
        param[i].result = &result[i];

        if(pthread_create(&th[i], NULL, thread, &param[i]) != 0){
            fprintf(stderr, "Failed to create thread%d.\n", i);
            exit(4);
        }

        printf("Thread%d has been created.\n", i);
    }

    //创建保存文件线程
    pthread_t saveth;
    unsigned int saveThreadFlag = 1;
    TParam saveParm;
    saveParm.cId = thread_num;      //保存线程序号紧跟在所有处理线程序号后
    saveParm.pOut = pOutBufArr;     //与处理线程不同，pout指向的不再是每个线程的起始位置，而是指向一个数组，数组中保存着各线程的缓冲区指针。
    saveParm.outSemEmpty = outSemEmpty;
    saveParm.outSemFull = outSemFull;
    saveParm.runFlag = &saveThreadFlag;

    if(pthread_create(&saveth, NULL, save, &saveParm) !=0){
        fprintf(stderr, "Failed to create save thread.\n");
        exit(4);
    }
    printf("Save thread has been created.\n");

     /*******************
    读取数据文件
    *******************/
    printf("%d\t%d\t%d\t%d\n", thread_num,fqNum,EACH_READ_NUM,READ_LEN);
    //    char  read_buf[thread_num][fqNum][EACH_READ_NUM][READ_LEN];         //不能用此方法申请内存，因为数组大小受限于栈的大小。默认程序的栈都有大小限制(据说编译器可调)

    //打开文件流的指针
    gzFile  RFILE[fqNum] ;
    for(int i=0; i<fqNum; i++){
        puts(fqA[i]);
        RFILE[i] = gzopen(fqA[i], "r");
        if(RFILE[i] == NULL){
            fprintf(stderr, "Error! Failed to open sequence data file: %s.\n", fqA[i]);
            exit(5);
        }
    }
    //循环读取文件内容
    char line[LINE_LEN];
    _Bool bFileEnd = false;
    int runTime =2 ;
    while(runTime){
        if(runTime==1){
            puts("+++++++++++ LAST +++++++++++");
        }
        for(int i=0; i<thread_num; i++){        //为每个线程读取数据

            //查看缓冲区是否为空
 //           printf("-------------MainThread is waiting (thread%d)%lu\n", i, &inSemEmpty[i]);
            sem_wait(&inSemEmpty[i]);
 //           printf("-------------ThreadMain has %lu\n", &inSemEmpty[i]);

            //最后一次循环时标记各处理线程的退出标志(不能把修改标记放到信号量处理外面)
            if(runTime==1){
                puts("最后一次空数据");
                printf("thread numer %d\n", i);
                //标记此线程退出标志
                pthread_mutex_lock(&mutex);
                t_flag[i] = 0;
                pthread_mutex_unlock(&mutex);
            }

            for(int j=0; j<fqNum&&(!bFileEnd); j++) {         //每端数据
                //memset(line, 0, LINE_LEN);            //
                for(int n=0; n<EACH_LINE_NUM&&(!bFileEnd); n++){     //读取指定行数的文件内容
                    if(gzgets(RFILE[j], line, LINE_LEN) == NULL){
                        //如果文件没有结束，则表示出错了
                        if(! gzeof(RFILE[j])){
                            fprintf(stderr, "Error: get null when read sequence data file(%s).\n", fqA[j]);
                            exit(6);
                        }
                        if(j == fqNum-1){       //所有fq文件已读完，但可能后面的几个线程没有分配到数据，所以还要另外单独处理
                            bFileEnd = true;

                            puts("File End!");
                        }else{                      //仅读完前面的fq文件
                            break;
                        }
                    }

                    memcpy(pInBufArr[i][j] +LINE_LEN*n, line, LINE_LEN);
                }
            }

            //while的判断变量减1
            if(bFileEnd && i==thread_num-1){
                puts("^^^^^^^ last time ^^^^^^");
                --runTime ;

            }

            //释放缓冲区已完成写入的信号
 //           printf("-------------ThreadMain will release %lu\n", &inSemFull[i]);
            sem_post(&inSemFull[i]);
 //           printf("-------------MainThread has released (thread%d)%lu\n", i, &inSemFull[i]);
        }
    }


    /*************
    等待子线程结束
    *************/
    for(int i = 0; i< thread_num; i++){
        pthread_join(th[i], NULL);
    }

    pthread_join(saveth, NULL);

    /**************
    主线程结束
    **************/
    puts("main thread finished.");
    return 0;
}

/***** 处理线程 *****/
void * thread(void *arg){
    TParam * pt = (TParam*)arg;
    printf("Thread %d is running...\n", pt->cId );
    unsigned  long buf_size = EACH_LINE_NUM*(LINE_LEN)*sizeof(char);

    //根据cid计算各信号量的指针位置
    sem_t  * pInSemEmpty = pt->inSemEmpty+ pt->cId;
    sem_t  * pInSemFull  = pt->inSemFull+ pt->cId;
    sem_t  * pOutSemEmpty = pt->outSemEmpty+ pt->cId;
    sem_t  * pOutSemFull = pt->outSemFull+ pt->cId;

 //   printf("########\nThread%u Sem(from thread):\ninEmpty: %lu\tinFull: %lu\toutEmpty: %lu\toutFull: %lu\n########\n", pt->cId, pInSemEmpty,pInSemFull,pOutSemEmpty,pOutSemFull);

    //缓冲存的指针
    char *pin[FQ_NUM];
    char *pout[FQ_NUM];

    //存放每条read的信息
    char readInfo[FQ_NUM][4][LINE_LEN];    //存放一条read信息（4行）
    memset(readInfo, 0, sizeof(readInfo));

    //统计结果存放变量
    TResult * result = pt->result;
    _Bool jumpOut = false;

    //循环等待接收数据并处理
    while(1){
//        printf("-------------Thread%u is waiting %lu\n", pt->cId, pOutSemEmpty);
        sem_wait(pOutSemEmpty);
//        printf("-------------ThreadMain has %lu\n", pOutSemEmpty);

 //       printf("-------------Thread%u is waiting %lu\n", pt->cId, pInSemFull);
        sem_wait(pInSemFull);
 //       printf("-------------Thread%u has %lu\n", pt->cId, pInSemFull);

        //将指针再指向缓冲区起始位置
        for(int i=0; i<FQ_NUM; i++){
            pin[i] = *(pt->pIn + i);
            pout[i] = *(pt->pOut + i);
        }

        //检查是否可退出(即是不是主线程最后一次空的数据)

//        puts("*******************");
//        printf("%u\t%d\n", *(pt->runFlag), *(*pt->pIn));
//        puts("*******************");
        if(*(pt->runFlag) == 0 && *(*pt->pIn) == '\0'){
            jumpOut = true;
            pthread_mutex_lock(&mutex);
            finished_thread_num++;
            pthread_mutex_unlock(&mutex);
            goto ignoreData;        //跳过后面的数据处理
        }


        //数据处理
        for(int i=0; i< EACH_READ_NUM; i++){    //reads
            //对读取的read进行分析。如果读取到的read为空，则表示缓冲区里的内容已经读完，跳出循环
            if(*(pin[0]) == '\0' || *(pin[1]) == '\0'){
                break;
            }

            result->allReadsNum += 1;

            for(int n=0; n< FQ_NUM; n++){           //单双端
                for(int j=0; j<4; j++){
                    memcpy(readInfo[n][j], pin[n], LINE_LEN);
                    pin[n] = pin[n]+LINE_LEN;   //移动指针到下一条read
             //       printf(readInfo[n][j]);
                }
            }
            //查看read情况
            //N值、低质量值、碱基含量
            _Bool bNDiscard = false;     //是否因为N较多去除read
            _Bool bLDiscard = false;      //是否因为低质量去除read

            int A[FQ_NUM*(READ_LEN-1)];
            int T[FQ_NUM*(READ_LEN-1)];
            int G[FQ_NUM*(READ_LEN-1)];
            int C[FQ_NUM*(READ_LEN-1)];
            int N[FQ_NUM*(READ_LEN-1)];

            int PHRD[FQ_NUM*(READ_LEN-1)];

            for(int x=0; x<FQ_NUM*(READ_LEN-1); x++){
                A[x] = 0;
                T[x] = 0;
                G[x] = 0;
                C[x] = 0;
                N[x] = 0;

                PHRD[x] = 0;
            }

 //           int phv[FQ_NUM][READ_LEN-1];        //read上每个位置的质量值

            for(int n =0; n<FQ_NUM ; n++){

                unsigned int phv = 0;
                unsigned int a_num = 0;
                unsigned int t_num = 0;
                unsigned int g_num = 0;
                unsigned int c_num = 0;
                unsigned int n_num = 0;

                for(int j=0; j< READ_LEN-1; j++){
                    int pos = n*(READ_LEN-1)+j;
                    //read
                    switch(readInfo[n][1][j] ){
                        case 'A':
                            A[pos] += 1;
                            a_num ++;
                            break;
                        case 'T':
                            T[pos] += 1;
                            t_num ++;
                            break;
                        case 'G':
                            G[pos] += 1;
                            g_num ++;
                            break;
                        case 'C':
                            C[pos] += 1;
                            c_num ++;
                            break;
                        case 'N':
                            N[pos] += 1;
                            n_num ++;
                            break;
                        default:
                            break;
                    }

                    phv += readInfo[n][3][j];
                    PHRD[pos] += readInfo[n][3][j];
                }
                //丢弃N较多的read
                if(n_num*1.0/(a_num +t_num+g_num+c_num) > N_Rate){
                    bNDiscard = true;
                    result->nCutoffReadsNum += 1;
                    break;
                }
                //-----------低质量------------
                if(phv*1.0/(READ_LEN-1) - PHRED_BASE < LOW_QUAL){
                    bLDiscard = true;
                    result->lowQualReadsNum += 1;
                    break;
                }
            }

           //如果是按整条read丢弃
           if(READ_MIN_BASES == 0) {
                if(bNDiscard || bLDiscard){
                    continue;       //忽略，读取下一条read
                }
                //统计clean data中各位置上各碱基数量和phred值
                for(int n=0; n<FQ_NUM*(READ_LEN-1); n++){
                    *(*(result->posBase + n) )  += A[n];
                    *(*(result->posBase + n) + 1)  += T[n];
                    *(*(result->posBase + n) + 2)  += G[n];
                    *(*(result->posBase + n) + 3)  += C[n];
                    *(*(result->posBase + n) + 4)  += N[n];

                    *(result->posPhred +n) += PHRD[n];
                }

                //把保留的reads数据写入相应的输出缓冲区

                for(int n=0; n<FQ_NUM; n++){
                    for(int j =0; j<4; j++){
                        memcpy(pout[n], readInfo[n][j], LINE_LEN);
                        pout[n] += LINE_LEN;
                 //       printf("%s",readInfo[n][j]);
                    }
                }

            //从3'端截短
           }else{

           }

         }

        //清空相应的缓冲区
       for(int i=0; i< FQ_NUM; i++){
           memset(*(pt->pIn + i), 0, buf_size);
       }


 //       printf("-------------Thread%u will release %lu\n", pt->cId, pInSemEmpty);
 ignoreData:       sem_post(pInSemEmpty);
 //       printf("-------------Thread%u  has released %lu\n", pt->cId, pInSemEmpty);
//        printf("-------------Thread%u will release %lu\n", pt->cId, pOutSemFull);
        sem_post(pOutSemFull);
 //       printf("-------------Thread%u has released %lu\n", pt->cId, pOutSemFull);

        if(jumpOut){
            break;
        }
    }

    //如果缓冲区数据不为空，说明runFlag被设置成0之前主线程读取的那次数据还未被本线程处理。
    printf("Thread%d finished!\n", pt->cId);
    return arg;
}

void * save(void *arg){
    unsigned  long buf_size = EACH_LINE_NUM*(LINE_LEN);        //每块缓冲区大小
    TParam * pt = (TParam*)arg;
    puts("Save thread is running...");
    const unsigned int FILE_NAME_LEN = NAME_LEN +64;
    char outfile[FQ_NUM][FILE_NAME_LEN];        //存放输出文件路径
    memset(outfile, 0, sizeof(outfile));
    for(int i=0; i<FQ_NUM; i++){
        snprintf(outfile[i], FILE_NAME_LEN, "%s/%s_%d.clean.fq",  out_dir, sample, i+1);
    }

    unsigned int thread_num = pt->cId;         //正好等于处理线程的数量
    _Bool brun = true;

    //打开文件流
    FILE * fileArr[FQ_NUM];
    for(int f=0; f<FQ_NUM; f++){
        fileArr[f] = fopen(outfile[f], "w");
    }

    //循环等待接收数据并写入最终的文件中
    while(brun){
        //按顺序检测信号量
        for(int tNum =0; tNum< thread_num; tNum++){
            sem_t  * pOutSemEmpty = pt->outSemEmpty +tNum  ;
            sem_t  * pOutSemFull = pt->outSemFull +tNum ;
   //         printf("==========\nThread%d Sem(from save):\noutEmpty: %lu\toutFull: %lu\n=========\n", tNum,pOutSemEmpty, pOutSemFull);

 //           printf("-------------SaveThread is waiting (thread%d)%lu\n", tNum, pOutSemFull);
            sem_wait(pOutSemFull);      //死锁？
 //           printf("-------------SaveThread has (thread%d)%lu\n", tNum, pOutSemFull);

            for(int n=0; n<FQ_NUM; n++){
                char * pstart = *(pt->pOut + FQ_NUM*tNum +n);
                char * pout = pstart ;      //指向处理线程对应的输出缓冲区
                //FILE * file = fopen(outfile[n], "w");
                //从缓冲区中读取数据并写入文件
                for(int i =0; i<EACH_READ_NUM; i++){
                    for(int j=0; j<4; j++){
                        fputs(pout, fileArr[n]);
                        pout += LINE_LEN;
                    }
                }
//                fflush(fileArr[n]);
                memset(pstart, 0, buf_size);
            }

            pthread_mutex_lock(&mutex);
            if(finished_thread_num == thread_num){
                brun = false;
            }
            pthread_mutex_unlock(&mutex);
//            printf("-------------Thread%u will release %lu\n", tNum, pOutSemEmpty);
            sem_post(pOutSemEmpty);
 //           printf("-------------SaveThread has released (thread%d)%lu\n", tNum, pOutSemEmpty);
        }

    }
    puts("SaveThread finished!");
    return arg;
}
