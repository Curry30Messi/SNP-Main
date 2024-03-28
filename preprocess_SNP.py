import xlrd

# what
def Read(file_path):
    book = xlrd.open_workbook(file_path)
    Table = book.sheet_by_index(0)
    return Table


class Output:
    def __init__(self, table,n):
        self.table = table
        self.normalChain = table.col_values(0)[1:210]
        self.antiChain = table.col_values(2)[1:210]
        self.name = table.col_values(3)[1:210]
        self.n = n
        self.dict1 = {}
        self.dict2 = {}
        self.dict3 = {}
        self.dict4 = {}

    def Search(self):
        pass

    def duplicateCheck(self, dict):  # 查重
        m = 0
        list = []
        for i in dict.values():
            list.append(i)
        list = sorted(list)
        print(list)
        for i in range(len(list)):
            for j in range(i + 1, len(list)):
                if list[j] > list[i]:
                    break
                else:
                    print(list[i], list[j])
                    m += 1
        return m

    def cutOut(self):  # 截取子串并编码
        for i in range(209):
            # 截取子串
            pos1 = self.normalChain[i].find('[') - 1
            pos2 = self.normalChain[i].find(']') + 2
            pos3 = self.antiChain[i].find('[') - 1
            pos4 = self.antiChain[i].find(']') + 2
            a = self.normalChain[i][pos1 - self.n:pos1]
            b = self.normalChain[i][pos2:pos2 + self.n]
            c = self.antiChain[i][pos3 - self.n:pos3]
            d = self.antiChain[i][pos4:pos4 + self.n + 1]
            # 编码并转化为十进制表示
            aa = self.binToDec(self.Encode(a))
            bb = self.binToDec(self.Encode(b))
            cc = self.binToDec(self.Encode(c))
            dd = self.binToDec(self.Encode(d))
            # 生成字典
            self.dict1[self.name[i] + ' + up'] = aa
            self.dict2[self.name[i] + ' + down'] = bb
            self.dict3[self.name[i] + ' - up'] = cc
            self.dict4[self.name[i] + ' - down'] = dd

    def Encode(self, str):  # 编码
        turn = ''
        for i in range(self.n):
            if str[i] == 'T':
                turn += '00'
            elif str[i] == 'A':
                turn += '01'
            elif str[i] == 'G':
                turn += '10'
            else:
                turn += '11'
        return turn

    def binToDec(self, bin):  # 二进制转十进制
        sum = 0
        for i in range(2 * self.n):
            sum += int(bin[2 * self.n - 1 - i]) * 2 ** i
        return sum



if __name__ == '__main__':
    table = Read('SNP/20230309--SNPs.xls')
    n = 8
    output = Output(table,n)
    output.cutOut()
    print(output.dict1)
    #print(output.duplicateCheck(output.dict2))
    #print(output.duplicateCheck(output.dict3))
    #print(output.duplicateCheck(output.dict4))



