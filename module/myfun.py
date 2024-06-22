import numpy as np
import pandas as pd
from functools import reduce

def saveExcel(A,filename):
    writer = pd.ExcelWriter(filename + '.xlsx')		# Creat Excel file
    if type(A) is list:
        for i,a in enumerate(A) :
            data = pd.DataFrame(a)   # convert to dataframe form
            data.to_excel(writer,'page '+str(i), float_format='%.5f')		# 写入excel,‘page_1’是写入excel的sheet名
    else:
        data = pd.DataFrame(A)   # convert to dataframe form
        data.to_excel(writer, float_format='%.5f')		# 写入excel,‘page_1’是写入excel的sheet名
    writer.save()  # Save file
    writer.close()  # Close file

# general code for getting web
def getHtmltxt(url):
    try:
        # get web
        r = requests.get(url)
        # 如果不是200，产生异常requests.HTTPError
        r.raise_for_status()
        #  从内容中分析出的响应内容编码方式
        r.encoding = r.apparent_encoding
        return r.text
    except:
        print("Error in getHtmlTxt")
        return ''

def functorial(n):
    '''Calculates n!'''
    return reduce( lambda x, y : x * y, (i+1 for i in range(n)) )

def permute(nums):
    ''' input: nums, a list of if elements
        return: full permutation of the input'''
    def permutation(nums, k, n): 
        if k == n:
            return res.append(nums[:])
        for i in range(k, n):  
            nums[i], nums[k] = nums[k], nums[i]  
            permutation(nums, k + 1, n) 
            nums[i], nums[k] = nums[k], nums[i] 
    res = []
    permutation(nums, 0, len(nums))
    return res
    
