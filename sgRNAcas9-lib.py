#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-11-21 20:15:20
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

'''
v1.0
'''

import os
import sys
from multiprocessing import Pool
import gc
import getopt
import platform
from bin import regex as re

def usage():
    print(
        """
Usage:
    python3 sgRNAcas9-lib.py -i <gene.fa> -g <genome.fa> -o <output.txt>
Options:
    -m    Maximum number of mismatches       <int>          [default:5]
    -p    The number of process              <int>          [default:8]
"""
    )


def getinput():
    global outputfile
    inputfile = ''
    outputfile = ''
    genomefile=''
    mismatch=5
    process=8
    if len(sys.argv) == 1:  # 无参数时 直接显示使用信息
        usage()
        sys.exit()
    else:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:g:p:m:", ["help"])  # 定义输入的- &--的情况
        for opt, value in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt == '-i':
                inputfile = value
                try:
                    open(inputfile)
                except:
                    print("please check your gene file")
                    sys.exit(2)
            elif opt == '-g':
                genomefile = value
                try:
                    open(genomefile)
                except:
                    print("please check your genome file")
                    sys.exit(2)
            elif opt == '-o':
                outputfile = value
                outputfile=outputfile.split('.')[0]+'.txt'
                if os.path.exists(outputfile):
                    os.remove(outputfile)
            elif opt == '-p':
                try:
                    process = int(value)
                except:
                    print("please check your process")
                    sys.exit(2)
            elif opt == '-m':
                try:
                    mismatch = int(value)
                except:
                    print("please check your mismatch")
                    sys.exit(2)
    with open(outputfile,'a+') as out:
        out.write('sgRNAid\tsgRNA\tstrand\tstart\tend\tGC\t')
        for i in range(mismatch+1):
            out.write('%sM\t'%(i))
        out.write('Total\tpostion\n')
    return inputfile,outputfile,genomefile,mismatch,process


def cleangenome(genomefile):
    #对参考基因组基因clean
    genomedict={}
    outfile=genomefile.split('.')[0]+'.clean.fa'
    if os.path.exists(outfile):
        pass
    else:
        out=open(outfile,'a+')
        with open(genomefile) as ff:
            for line in ff:
                line=line.strip()
                if line.startswith('>'):
                    name=line.split('.')[0].replace('>','')
                    genomedict[name]=[]
                else:
                    genomedict[name].append(line.upper())
        for key,value in genomedict.items():
                out.write('>'+key+'\n')
                value=''.join(value)
                out.write(value+'\n')
        out.close()
    return outfile


def CountGCcontent(seq):
    seq = seq.upper() #也可是使用.lower()把大写转换成小写计算
    count_c = seq.count('C')
    count_g = seq.count('G')
    gc_content = ((count_g + count_c) / len(seq))*100
    return round(gc_content,2)


def calTTTT(seq):
    seq=seq.upper()
    num1=seq.count('TTTT')
    num2=len(re.findall(r'([A-Z])\1{5}',seq))
    num=num1+num2
    if num>0:
        return False
    else:
        return True

def Fasta_reverse(sequence):
    #将序列进行方向互补
    sequence=sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.upper()
    return sequence[::-1]


def output(offList):
    with open(outputfile,'a+') as out:
        outinfo='\t'.join(offList)
        out.write(outinfo+'\n')


def calOffTarget(Results,name,seq,strand,start,end,mismatch):
    os.system('%s -i %s -p NGG -n %s -g mm10 -o %s'%(Results,seq,mismatch,name))
    gcgc=CountGCcontent(seq[:20])
    offnumlist=[0 for i in range(mismatch+1)]
    offList=[name,seq,strand,start,end,gcgc]
    postion=[]
    with open(name) as ff:
        for line in ff:
            line=line.strip()
            offnum=int(line.split()[-1])
            if offnum==0:
                test=line.split()[:4]
                postion.append('_'.join(test))
            offnumlist[offnum]+=1
        total=sum(offnumlist)-1
        if total<0:
            total=0
        offnumlist.append(total)
        if len(postion)>0 and len(postion)<10:
            postion='|'.join(postion)
            offnumlist.append(postion)
        elif len(postion)>=10:
            postion='too many off-target'
            offnumlist.append(postion)
        else:
            postion='error'
            offnumlist.append(postion)
    offList.extend(offnumlist)
    offList=[str(i) for i in offList]
    os.remove(name)
    return offList
    del gcgc,postion,offnum,offList
    gc.collect()

def run(genefile,Results,mismatch=5,process=8):
    geneDict={}
    with open(genefile) as ff:
        for line in ff:
            line=line.strip()
            if line.startswith('>'):
                name=line.replace('>','')
                geneDict[name]=[]
            else:
                line=line.upper()
                geneDict[name].append(line)
    p=Pool(process)
    pattern1 = re.compile('.{21}GG')
    pattern2 = re.compile('CC.{21}')
    for key in list(geneDict.keys()):
        value=''.join(geneDict[key])
        del geneDict[key]
        sgRNA_L=pattern1.finditer(value,overlapped=True)
        sgRNA_F=pattern2.finditer(value,overlapped=True)
        count=0
        for sgRNA in sgRNA_L:
            count+=1
            name=key+'_s_'+str(count)
            start=sgRNA.start()
            end=sgRNA.end()
            seq=sgRNA.group()
            if calTTTT(seq[:20]):
                p.apply_async(calOffTarget,args=(Results,name,seq,'+',start,end,mismatch),callback=output)
        count=0
        for sgRNA in sgRNA_F:
            count+=1
            name=key+'_a_'+str(count)
            start=sgRNA.start()
            end=sgRNA.end()
            seq=Fasta_reverse(sgRNA.group())
            if calTTTT(seq[:20]):
                p.apply_async(calOffTarget,args=(Results,name,seq,'-',start,end,mismatch),callback=output)
    p.close()
    p.join()

def del_files(path_file):
    ls = os.listdir(path_file)
    for i in ls:
        f_path = os.path.join(path_file, i)
        # 判断是否是一个目录,若是,则递归删除
        if os.path.isdir(f_path):
            del_files(f_path)
        else:
            os.remove(f_path)

def check_bin():
    if os.path.exists('bin'):
        pass
    else:
        print('The bin file is missing')
        sys.exit(2)


def check_install(Table_Creation):
    if os.path.exists(Table_Creation):
        pass
    else:
        os.system('cd bin&&make')

def check_index(Table_Creation,genomefile,genomedirr):
    if os.path.exists(genomedirr):
        if os.path.exists(genomedirr+genomefile):
            print('index file already exists')
        else:
            del_files(genomedirr)
            os.system('%s -i %s -o %s'%(Table_Creation,genomefile,genomedirr))
    else:
        os.makedirs(genomedirr)
        os.system('%s -i %s -o %s'%(Table_Creation,genomefile,genomedirr))



def main(inputfile,outputfile,genomefile,mismatch,process):
    print('Initializing file,please wait!')
    genomefile=cleangenome(genomefile)
    ##################
    dirr='./bin/'
    genomedirr='./genome/'
    ##################
    Table_Creation=dirr+'Table_Creation'
    Results=dirr+'Results'
    Load_Memory=dirr+'Load_Memory'
    Detach_Memory=dirr+'Detach_Memory'
    ##################
    check_install(Table_Creation)
    check_index(Table_Creation,genomefile,genomedirr)
    ######导入内存######
    os.system('%s -t %s -g mm10'%(Load_Memory,genomedirr))
    #####结果##########
    run(inputfile,Results,mismatch=mismatch,process=process)
    print('Detaching Memory,please wait!')
    ##################
    os.system('%s -g mm10'%(Detach_Memory))
    ##################
    out=open(genomedirr+genomefile,'a+')
    out.close()



if __name__ == '__main__':
    if(platform.system()=='Linux'):
        check_bin()
        inputfile,outputfile,genomefile,mismatch,process=getinput()
        main(inputfile,outputfile,genomefile,mismatch,process)
    else:
        print('Please use linux system to run the program')
        sys.exit(2)
