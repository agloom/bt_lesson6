#ifndef QC_H_INCLUDED
#define QC_H_INCLUDED

void * thread(void *arg);
void * save(void *arg);
/**
根据fq内容判断phred33还是phred64
*/
int phredVal(char * fq_file_path, unsigned int lineLen);

typedef struct _thread_result{
    unsigned long long aNum ;
    unsigned long long tNum;
    unsigned long long gNum ;
    unsigned long long cNum ;
    unsigned long long nNum ;

    unsigned long long lowQualReadsNum ;     //低质量reads数
    unsigned long long nCutoffReadsNum ;     //N大于cutoff的reads数
    unsigned long long allReadsNum ;             //处理的reads总数

    unsigned long long *posPhred;       //指向一个数组的起始位置，每个数组元素表示此位置上碱基质量值的和
    unsigned long long (*posBase)[5];        //指向包含5个元素数据的longlong指针，5个元素分别表示atgcn碱基
    unsigned long long q20Num ;                      //Q20以上的碱基数
    unsigned long long q30Num ;                      //Q30以上的碱基数
} TResult;

typedef struct _TParam{
    unsigned int cId;       //创建的第几个子线程，从0开始
    char ** pIn;                 //存放要读取数据的缓冲区起始位置
    char ** pOut;               //存放要存入数据
    sem_t * inSemEmpty;   //读入数据的锁
    sem_t * inSemFull;
    sem_t * outSemEmpty;   //写出数据的锁
    sem_t * outSemFull;
    unsigned int * runFlag;     //是否保持运行，1表示保持运行，0表示可退出。
    TResult * result;                  //存放线程运行结果的数据结构
} TParam;

#endif // QC_H_INCLUDED
