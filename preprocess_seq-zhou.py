import pandas as pd
import os


def Chain(bef):
    if bef[-2:] == '/2':
        return '/2'
    elif bef[-2:] == '/1':
        return '/1'
    else:
        return '00'


def seq_standardize():
    files = os.listdir('data')
    for filename in files:
        dict = {}
        ext = os.path.splitext(filename)[1]     # 取文件后缀名
        if ext == ".fastq":
            with open('data/' + filename) as fq:
                bef = ''
                for line in fq:
                    # one
                    br = 1
                    line = line.replace('\n', '')
                    if line == '':
                        continue
                    for j in range(len(line)):
                        if line[j] in 'ACGTN':
                            pass
                        else:
                            br = 0
                            break
                    if br:
                        line0 = bef[1:].split(' ')
                        seq_name = line0[0]
                        dict[seq_name] = line
                        chain = Chain(bef)
                        dict[seq_name] += chain
                    bef = line
        elif ext == ".fa" or ext == ".fna":
            # print(1)
            with open('data/' + filename) as fq:
                fq = fq.read()
                number = 1
                first = fq.find('>')
                while True:
                    start = fq.find('\n', first)
                    end = fq.find('>', start)
                    seq_name = str(number)
                    number += 1
                    dict[seq_name] = fq[start:end].replace('\n', '') + "/1"
                    if end != -1:
                        first = end
                    else:
                        break
                # # two
                # if line == '\n':
                #     # 处理空白行情况
                #     i = 0
                #     continue
                # if i % 4 == 1:
                #     # 去除末尾换行符
                #     line = line.replace('\n', '')
                #     line = line[1:]
                #     line0 = line.split(' ')
                #     seq_name = line0[0]
                #     chain = line[-1]
                #     # print(seq_name)
                #     dict[seq_name] = ''
                # elif i % 4 == 2:
                #     seq = line.replace('\n', '')
                #     dict[seq_name] = seq
                # elif i % 4 == 3:
                #     line = line.replace('\n', '')
                #     feature = line
                #     dict[seq_name] += chain
        # 将结果写入csv文件中
        result = pd.Series(dict)
        # result.name = feature
        result.index.name = 'ip'
        result.to_csv('data_std/' + filename[:-6] + '.csv', index=True, encoding='gbk')
    print("End of standardization\n")


if __name__ == '__main__':
    seq_standardize()
