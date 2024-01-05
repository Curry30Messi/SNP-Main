from tqdm import tqdm
import csv
import openpyxl
import xlrd
import os
import time
from preprocess_seq import seq_standardize
import random
from openpyxl import Workbook
from openpyxl import load_workbook
import pandas as pd
from openpyxl.utils import get_column_letter
from collections import Counter


def countFile(output, file_path):
    with open(file_path, 'w', newline='') as csvfile:
        w = csv.writer(csvfile)
        w.writerow(['SNPname', 'A', 'T', 'G', 'C'])
        if output.feature == '/1':
            for i in range(len(output.normal)):
                w.writerow([output.normal[i][0], output.normal[i][3], output.normal[i][4], output.normal[i][5], output.normal[i][6]])
        else:
            for i in range(len(output.anti)):
                w.writerow([output.anti[i][0], output.anti[i][3], output.anti[i][4], output.anti[i][5], output.anti[i][6]])


def Read(file_path):
    book = xlrd.open_workbook(file_path)
    Table = book.sheet_by_index(0)
    return Table


def KMP(mon_string, son_string):
    i = 0
    j = 0
    next = get_next(son_string)
    while (i < len(mon_string) and j < len(son_string)):
        if j == -1 or mon_string[i] == son_string[j]:
            i += 1
            j += 1
        else:
            j = next[j]

    if j == len(son_string):
        return i - j
    return -1


def get_next(son_string):
    next = [-1] * len(son_string)
    next[1] = 0
    i = 1
    j = 0
    while i < len(son_string) - 1:
        if j == -1 or son_string[i] == son_string[j]:
            i += 1
            j += 1
            next[i] = j
        else:
            j = next[j]

    return next


def Anti(vlu):
    if vlu == 'T':
        vlu = 'A'
    elif vlu == 'A':
        vlu = 'T'
    elif vlu == 'G':
        vlu = 'C'
    elif vlu == 'C':
        vlu = 'G'
    return vlu
# accu
def seq_read(output, result_file, time1, snp_count_stand, q):
# fuzzy
# def seq_read(output, result_file, time1, cnumber):
    n = output.max_k  # k为上下游单侧范围取值
    files = os.listdir('data_std')
    with open(result_file, 'w', newline='') as csvfile:
        w = csv.writer(csvfile)
        # accu
        w.writerow(['id', 'SNP——IP', '突变位点碱基', '突变位置', '单侧匹配长度', '两侧匹配长度', '正反链'])
        # fuzzy
        # w.writerow(['id', 'SNP——IP', '单侧上游匹配度', '单侧下游匹配度', '两侧上游匹配度', '两侧下游匹配度', '突变位点碱基', '突变位置', '单侧匹配长度', '两侧匹配长度', '正反链'])

    for filename in files:
        ser = pd.read_csv('data_std/' + filename)
        # ser = ser.sample(frac=1.0)  # 随机打乱所有数据
        clo_list = ser.columns
        whole_gene = pd.Series(ser[clo_list[1]].values, index=ser[clo_list[0]].values)

    # print(len(whole_gene), ' Total')
    snp_total = output.normal + output.anti
    i = 0   # 用于记录SNP的
    for SNP in tqdm(snp_total):
        snp_count = 0
        time_start = time.time()
        ok = 0
        result_item = []
        for id, seq in whole_gene.items():
            # fuzzy
            if time.time() - time_start > time1 or snp_count >= snp_count_stand:
                break
            output.feature = seq[-2:]
            seq = seq[:-2]
            # fuzzy
            # if output.oneortwo == 1:
            #     if output.feature == '/1' or output.feature == '00':
            #         match = output.normalMatch1
            #     elif output.feature == '/2':
            #         match = output.antiMatch1
            # if output.oneortwo == 2:
            #     if output.feature == '/1' or output.feature == '00':
            #         match = output.normalMatch2
            #     elif output.feature == '/2':
            #         match = output.antiMatch2
            # result = []
            # id:全基因组中的片段编号
            # seq:id所对应的片段序列，待处理
            ln = len(seq)
            # KMPq
            list = []
            q_loc = KMP(seq, SNP[1][-q:])
            if q_loc != -1 and q_loc - len(SNP[1]) + q >= 0:
                if seq[q_loc - len(SNP[1]) + q:q_loc + q] == SNP[1] and len(seq) >= q_loc + q + 1:
                    snp_count += 1
                    list.append(id)
                    list.append(SNP[0])
                    if output.feature == '/1' or output.feature == '00':
                        list.append(seq[q_loc + q])
                        output.Count(i % 239, seq[q_loc + q])
                    elif output.feature == '/2':
                        list.append(Anti(seq[q_loc + q]))
                        output.Count(i % 239, Anti(seq[q_loc + q]))
                    list.append(str(q_loc + q + 1))
                    list.append(str(output.k1))
                    list.append(str(output.k2))
                    list.append(output.feature)
                    with open(result_file, 'a+', newline='') as csvfile:
                        w = csv.writer(csvfile)
                        w.writerow(list)
                # a = KMP(seq[q_loc - len(SNP[1]) + q:], SNP[1])
                # b = KMP(seq, SNP[2])
                # # c = KMP(seq, SNP[1][-output.k2:])
                # # d = KMP(seq, SNP[2][:output.k2])
                # if a != -1 and a + output.k1 + 1 <= len(seq):
                #     snp_count += 1
                #     list.append(id)
                #     list.append(SNP[0])
                #     if output.feature == '/1' or output.feature == '00':
                #         list.append(seq[a + output.k1])
                #     elif output.feature == '/2':
                #         list.append(Anti(seq[a + output.k1]))
                #     list.append(str(a + output.k1 + 1))
                #     list.append(str(output.k1))
                #     list.append(str(output.k2))
                #     list.append(output.feature)
                #     with open(result_file, 'a+', newline='') as csvfile:
                #         w = csv.writer(csvfile)
                #         w.writerow(list)
                # if b != -1:
                #     snp_count += 1
                #     list.append(id)
                #     list.append(SNP[0])
                #     if output.feature == '/1' or output.feature == '00':
                #         list.append(seq[b - 1])
                #     elif output.feature == '/2':
                #         list.append(Anti(seq[b - 1]))
                #     list.append(str(b))
                #     list.append(str(output.k1))
                #     list.append(str(output.k2))
                #     list.append(output.feature)
                #     with open(result_file, 'a+', newline='') as csvfile:
                #         w = csv.writer(csvfile)
                #         w.writerow(list)

        i += 1
    # process_table2(filename)
        # time.sleep(0.05)

class Output:
    def __init__(self, table, k1, k2):
    # def __init__(self, table, k1, k2, unilateral, bilateral, oneortwo):    # 传入SNP表格，k值，单侧匹配度要求，两侧匹配度要求
        self.table = table
        self.normalChain = table.col_values(2)[1:240]
        self.antiChain = table.col_values(3)[1:240]
        self.name = table.col_values(1)[1:240]
        self.k1 = k1
        self.k2 = k2
        self.max_k = max(k1, k2)
        self.feature = ''
        self.normal = []
        self.anti = []
        # self.unilateral = unilateral    # 单侧匹配度要求
        # self.bilateral = bilateral      # 两侧匹配度要求
        # self.time = time
        # self.oneortwo = oneortwo

    def cutOut(self):  # 截取子串并存储
        for i in range(239):
            list1 = []
            list2 = []
            # 截取子串
            pos1 = self.normalChain[i].find('[')
            pos2 = self.normalChain[i].find(']') + 1
            pos3 = self.antiChain[i].find('[')
            pos4 = self.antiChain[i].find(']') + 1
            a = self.normalChain[i][pos1 - self.max_k:pos1]
            b = self.normalChain[i][pos2:pos2 + self.max_k]
            c = self.antiChain[i][pos3 - self.max_k:pos3]
            d = self.antiChain[i][pos4:pos4 + self.max_k]
            # 将正链和反链的IP，上游和下游分别存入self.normal和self.anti两个list中
            list1.append(self.name[i])
            list1.append(a)
            list1.append(b)
            list1.append(0)     # 统计碱基A个数
            list1.append(0)     # 统计碱基T个数
            list1.append(0)     # 统计碱基G个数
            list1.append(0)     # 统计碱基C个数
            self.normal.append(list1)
            list2.append(self.name[i])
            list2.append(c)
            list2.append(d)
            list2.append(0)     # 统计碱基A个数
            list2.append(0)     # 统计碱基T个数
            list2.append(0)     # 统计碱基G个数
            list2.append(0)     # 统计碱基C个数
            self.anti.append(list2)

    def Count(self, i, vlu):
            if vlu == 'A':
                self.normal[i][3] += 1
            if vlu == 'T':
                self.normal[i][4] += 1
            if vlu == 'G':
                self.normal[i][5] += 1
            if vlu == 'C':
                self.normal[i][6] += 1

    def normalMatch1(self, str1, str2, id, SNP, vlu, pos):    # 传入匹配序列的上下游子串，全基因组的片段序号，SNP, 突变点碱基，突变点位置(正链用匹配函数)
        match = []  # 存储符合匹配度要求的IP
        list1 = []
        s1_up = self.SimilarityAccu(SNP[1][-k1:], str1[-k1:])
        s1_down = self.SimilarityAccu(SNP[2][:k1], str2[:k1])
        s2_up = self.SimilarityAccu(SNP[1][-k2:], str1[-k2:])
        s2_down = self.SimilarityAccu(SNP[2][:k2], str2[:k2])
        k = 0
        if s1_up or s1_down or (s2_up and s2_down):
            # self.Count(i, vlu)
            list1.append(id)
            list1.append(SNP[0])
            list1.append(vlu)
            list1.append(pos)
            list1.append(str(self.k1))
            list1.append(str(self.k2))
            list1.append(self.feature)
            match.append(list1)
            k = 1
        # s1_up = self.SimilarityFuzzy(SNP[1][-k1:], str1[-k1:], k1)
        # s1_down = self.SimilarityFuzzy(SNP[2][:k1], str2[:k1], k1)
        # k = 0
        # if s1_up >= self.unilateral or s1_down >= self.unilateral:
        #     # self.Count(i, vlu)
        #     list1.append(id)
        #     list1.append(SNP[0])
        #     list1.append(str(s1_up))
        #     list1.append(str(s1_down))
        #     list1.append(str(''))
        #     list1.append(str(''))
        #     list1.append(vlu)
        #     list1.append(pos)
        #     list1.append(str(self.k1))
        #     list1.append(str(self.k2))
        #     list1.append(self.feature)
        #     match.append(list1)
        #     k = 1
        return match, k    # 返回符合匹配度的SNP的IP以及上下游匹配度

    def normalMatch2(self, str1, str2, id, SNP, vlu, pos):    # 传入匹配序列的上下游子串，全基因组的片段序号，SNP, 突变点碱基，突变点位置(正链用匹配函数)
        match = []  # 存储符合匹配度要求的IP
        list1 = []
        s2_up = self.SimilarityFuzzy(SNP[1][-k2:], str1[-k2:], k2)
        s2_down = self.SimilarityFuzzy(SNP[2][:k2], str2[:k2], k2)
        # if s1_up > 0.6:
        #     print(s1_up)
        k = 0
        if s2_up + s2_down >= 2 * self.bilateral:
            # self.Count(i, vlu)
            list1.append(id)
            list1.append(SNP[0])
            list1.append(str(''))
            list1.append(str(''))
            list1.append(str(s2_up))
            list1.append(str(s2_down))
            list1.append(vlu)
            list1.append(pos)
            list1.append(str(self.k1))
            list1.append(str(self.k2))
            list1.append(self.feature)
            match.append(list1)
            k = 1
        return match, k    # 返回符合匹配度的SNP的IP以及上下游匹配度

    def antiMatch1(self, str1, str2, id, SNP, vlu, pos):    # 传入匹配序列的上下游子串    (反链用匹配函数)
        match = []  # 存储符合匹配度要求的IP
        if vlu == 'T':
            vlu = 'A'
        elif vlu == 'A':
            vlu = 'T'
        elif vlu == 'G':
            vlu = 'C'
        elif vlu == 'C':
            vlu = 'G'
        list1 = []
        s1_up = self.SimilarityAccu(SNP[1][-k1:], str1[-k1:])
        s1_down = self.SimilarityAccu(SNP[2][:k1], str2[:k1])
        s2_up = self.SimilarityAccu(SNP[1][-k2:], str1[-k2:])
        s2_down = self.SimilarityAccu(SNP[2][:k2], str2[:k2])
        k = 0
        if s1_up or s1_down or (s2_up and s2_down):
            # self.Count(i, vlu)
            list1.append(id)
            list1.append(SNP[0])
            list1.append(vlu)
            list1.append(pos)
            list1.append(str(self.k1))
            list1.append(str(self.k2))
            list1.append(self.feature)
            match.append(list1)
            k = 1
        return match, k  # 返回符合匹配度的SNP的IP以及上下游匹配度

    def antiMatch2(self, str1, str2, id, SNP, vlu, pos):    # 传入匹配序列的上下游子串    (反链用匹配函数)
        match = []  # 存储符合匹配度要求的IP
        if vlu == 'T':
            vlu = 'A'
        elif vlu == 'A':
            vlu = 'T'
        elif vlu == 'G':
            vlu = 'C'
        elif vlu == 'C':
            vlu = 'G'
        list1 = []
        s2_up = self.SimilarityFuzzy(SNP[1][-k2:], str1[-k2:], k2)
        s2_down = self.SimilarityFuzzy(SNP[2][:k2], str2[:k2], k2)
        k = 0
        if s2_up + s2_down >= 2 * self.bilateral:
            # self.Count(i, vlu)
            list1.append(id)
            list1.append(SNP[0])
            list1.append(str(''))    # 上下游匹配度交换
            list1.append(str(''))
            list1.append(str(s2_down))
            list1.append(str(s2_up))
            list1.append(vlu)
            list1.append(pos)
            list1.append(str(self.k1))
            list1.append(str(self.k2))
            list1.append(self.feature)
            match.append(list1)
            k = 1
        return match, k  # 返回符合匹配度的SNP的IP以及上下游匹配度

    def SimilarityFuzzy(self, str1, str2, k):   # 传入要进行相似度计算的两子串
        count = 0
        if self.unilateral == 1:
            if str1 == str2:
                return 1
        else:
            for i in range(k):
                if str1[i] == str2[i]:
                    count += 1
        return count / k   # 返回相似度

    def SimilarityAccu(self, str1, str2):  # 传入要进行相似度计算的两子串
        return str1 == str2


if __name__ == '__main__':
    result_file = 'result.csv'
    count_file = 'count.csv'
    table = Read('SNP/NewSNP.xls')
    f = open('requirement_fuzzy.txt')
    # accu
    k1 = int(f.readline().replace('\n', ''))
    # k2 = int(f.readline().replace('\n', ''))
    time1 = float(f.readline().replace('\n', ''))
    snp_count_stand = int(f.readline().replace('\n', ''))
    q = int(f.readline().replace('\n', ''))
    k2 = 0
    print('单侧匹配长度：', k1, '一条snp搜索时间:', time1, '一条snp的搜索次数:', snp_count_stand, '预匹配长度q:', q)
    output = Output(table, k1, k2)
    # fuzzy
    # k1 = int(f.readline().replace('\n', ''))
    # k2 = int(f.readline().replace('\n', ''))
    # ul = float(f.readline().replace('\n', ''))
    # bl = float(f.readline().replace('\n', ''))
    # time1 = float(f.readline().replace('\n', ''))
    # cnumber = int(f.readline().replace('\n', ''))
    # oneortwo = int(f.readline().replace('\n', ''))
    # print('单侧匹配长度：', k1, '两侧匹配长度：', k2,  '单侧匹配度要求：', ul, '两侧匹配度要求：', bl, '截至时长：', time1, 'hours',
    #       '检测次数：', cnumber, '单双侧检测方式：', oneortwo)
    # output = Output(table, k1, k2, ul, bl, oneortwo)
    output.cutOut()
    print('end of cutout')
    seq_standardize()
    # accu
    seq_read(output, result_file, time1, snp_count_stand, q)
    # fuzzy
    # seq_read(output, result_file, time1 * 3600, cnumber)
    countFile(output, count_file)
    print('end of matchingpip')