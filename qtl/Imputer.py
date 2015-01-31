'''
Created on Dec 4, 2014

@author: jiao
'''

def impute(list,**kwargs):
    '''
    It imputes the missing value in an array.
    '''
    missing_values = 'NaN' # default missing value type.
    if kwargs.get('missing_values'): 
        missing_values = kwargs.get('missing_values')
    for i in range(len(list)):
        prev = None
        next = None
        gap = None
        step = None
        if list[i] == missing_values:
            for j in range(i,len(list)):     
                if list[j] != missing_values:
                    if i!=0:
                        prev = i-1
                        next = j
                        gap = j-i
                        #print 'next',list[next],'prev',list[prev]
                        step = (list[next]- list[prev])/(gap+1.0)
                        #print step, gap,'not end'
                        for k in range(gap):
                            list[i+k] = list[i-1+k]+step
                        break
                    else:
                        list[0] = 0
                        prev = i
                        next = j
                        gap = j-i
                        step = (list[next]- list[prev])/float(gap)
                        for k in range(1,gap):
                            list[i+k] = list[i-1+k]+step
                        break                  
                else:
                    if j ==len(list) -1:
                        list[j] = 0
                        prev = i-1
                        gap = j-i
                        step = (list[j]-list[prev])/(gap+1.0)
                        #print gap,step,'end'
                        for k in range(gap-1):
                            list[i+k] = list[i-1+k]+step
                        break
                    else:
                        pass
        else:
            pass           
    return list

if __name__ == '__main__':
    pass
    #list = [1,2,3,4,5,6,'NaN','NaN',9]
    #list1 = [1,2,3,4,5,6,'NaN','NaN','NaN']
    #list2 = ['NaN','NaN','NaN',3,4,5,6,7,8,9]
    #list3 = ['NaN','NaN',1,2,3,4,5,6,7,8,9,10,'NaN','NaN',13,14,15,16,17,18,'NaN','NaN','NaN']
    #impute(list3,missing_values = 'NaN')