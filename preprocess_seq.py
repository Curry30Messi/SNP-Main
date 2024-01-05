import pandas as pd
import os


def Chain(bef):
    if bef[-2:] == '/2':
        return '/2'
    elif bef[-2:] == '/1':
        return '/1'
    else:
        return '00'


def check_sequence(input_str):
    valid_chars = set('ATCGN')
    chars = set(input_str)

    if chars.issubset(valid_chars) and len(chars.intersection('ATCG')) == 4:
        return True
    else:
        return False


def check_sequence_last(input_str):
    valid_chars = set('atcgnATCGN0123456789/')
    chars = set(input_str)

    if chars.issubset(valid_chars) and len(chars.intersection('ATCG')) == 4:
        return True
    else:
        return False


def seq_standardize(folder):
    files = os.listdir(folder)
    for filename in files:
        dict = {}
        current_key = None
        ext = os.path.splitext(filename)[1]  # 取文件后缀名
        if ext == ".fastq":
            with open(folder + filename) as fq:
                num_error = 1
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
                        if "." in seq_name and any(char.isdigit() for char in seq_name):
                            # 如果seq_name包含句点且至少包含一个数字，则不进行任何操作
                            pass
                        else:
                            # 如果seq_name不符合条件，将其重新赋值为错误消息
                            temp1=num_error
                            seq_name = "not_standard_file_error___" + str(temp1)
                            num_error += 1

                        dict[seq_name] = line
                        chain = Chain(bef)
                        dict[seq_name] += chain

                    bef = line
        elif ext == ".fa" or ext == ".fna":
            # print(1)
            with open(folder + filename) as fq:
                fq = fq.read()
                number = 1
                first = fq.find('>')
                while True:
                    start = fq.find('\n', first)
                    end = fq.find('>', start)
                    seq_name = str(number)
                    number += 1
                    dict[seq_name] = fq[start:end].replace('\n', '') + "/0"
                    if end != -1:
                        first = end
                    else:
                        break

        # 将结果写入csv文件中
        # my_dict = {key: value for key, value in dict.items() if check_sequence_last(value)}
        to_delete = []  # 创建一个临时列表，用于存储需要删除的键

        for key, value in dict.items():
            if check_sequence_last(value):
                pass
            else:
                to_delete.append(key)  # 将需要删除的键添加到临时列表中

        # 在循环结束后根据临时列表删除字典中的键值对
        for key in to_delete:
            dict.pop(key)
        result = pd.Series(dict)
        # result.name = feature
        result.index.name = 'ip'
        if filename[-6:] == ".fastq":
            result.to_csv('data_std/' + filename[:-6] + '.csv', index=True, encoding='gbk')
        elif filename[-3:] == ".fa":
            result.to_csv('data_std/' + filename[:-3] + '.csv', index=True, encoding='gbk')
        elif filename[-4:] == ".fna":
            result.to_csv('data_std/' + filename[:-4] + '.csv', index=True, encoding='gbk')
    print("End of standardization\n")


if __name__ == '__main__':
    folder_name = input("输入文件夹名称地址")
    seq_standardize(folder_name)
