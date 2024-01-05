import os
import csv
import pandas as pd
from pandas import Series
from preprocess_SNP import Output, Read
from preprocess_seq import seq_standardize

def seq_read(output,result_file):
    n = output.n  # k为上下游单侧范围取值
    files = os.listdir('data_std')
    with open(result_file, 'w',newline='') as csvfile:
        w = csv.writer(csvfile)
        w.writerow(['id','sequence','mutated nucleotide','location','match range(oneside)'])

    for filename in files:
        d_for = {}
        d_res = {}
        ser = pd.read_csv('data_std/' + filename)
        clo_list = ser.columns
        # print(clo_list)
        # print(type(ser[clo_list[1]]))
        whole_gene = pd.Series(ser[clo_list[1]].values,index=ser[clo_list[0]].values)
        # print(whole_gene.items)

        feature = clo_list[1]  # feature('+'/'-'):全基因组片段为正链或反

        if feature == '+':
            d_for = output.dict1
            d_res = output.dict2
        elif feature == '-':
            d_for = output.dict3
            d_res = output.dict4
        # print("dict")
        d_for_keys = list(d_for.keys())
        d_for_vals = list(d_for.values())
        d_res_keys = list(d_res.keys())
        d_res_vals = list(d_res.values())
        for id, seq in whole_gene.items():
            result = []
            # id:全基因组中的片段编号
            # seq:id所对应的片段序列，待处理
            ln = len(seq)
            # 取正链
            number = 0
            for i in range(0, ln - 2*n):
                number += 1
                # print('i:'+str(i))
                seq_for = seq[i:i + n]
                # 取正链(forward strand)长度为n,反链(reverse strand)
                seq_for_bin = output.Encode(seq_for)
                seq_for_Dec = output.binToDec(seq_for_bin)
                # 在正链突变位点字典中寻找并返回
                flag_0 = 1
                if seq_for_Dec in d_for_vals:
                    ind1 = d_for_vals.index(seq_for_Dec)
                    key = d_for_keys[ind1]
                    flag_0 = 0
                # for name_for,value_for in d_for.items():
                #     if seq_for_Dec == value_for:
                #         key = name_for
                #         flag_0 = 0
                if flag_0:
                    # print('no find*'+str(number))
                    continue

                # print('after key\n')
                # key = 'name + up'
                key_list = key.split(feature)
                key_rev = key_list[0] + feature +' down'
                flag_1 = 1
                ind2 = d_res_keys.index(key_rev)
                value_rev = d_res_vals[ind2]
                flag_1 = 0
                # for name_res,value_res in d_res.items():
                #     if key_rev == name_res:
                #         value_rev = value_res
                #         flag_1 = 0
                #         break
                if flag_1:
                    continue
                seq_rev = seq[i+n+1:i+2*n+1]

                seq_rev_bin = output.Encode(seq_rev)
                seq_rev_Dec = output.binToDec(seq_rev_bin)

                if value_rev == seq_rev_Dec:
                    # 匹配成功
                    result = [id,key_list[0],seq[i+n],str(i+n+1),str(n)]
                    # 全基因组片段突变点位置信息str(i+n+1)从1开始
                    # print('1')
                    with open(result_file,'a+',newline = '') as csvfile:
                        w = csv.writer(csvfile)
                        w.writerow(result)
                    break

if __name__ == '__main__':
    result_file = 'result_accu.csv' # 匹配结果文件位置
    SNP_file = 'SNP/20230309--SNPs.xls' # SNP文件位置
    f = open('requirement_accu.txt')
    n = int(f.readline())

    table = Read(SNP_file)
    output = Output(table,n)
    output.cutOut()
    print('end of cutout')
    # seq_standardize()
    seq_read(output,result_file)

    print('end of matching')