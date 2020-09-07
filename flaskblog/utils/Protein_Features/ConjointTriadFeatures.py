#!/usr/bin/env python
# coding: utf-8

# In[1]:


import string


# In[2]:


AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]


# In[3]:


_repmat={1:["A",'G','V'],2:['I','L','F','P'],3:['Y','M','T','S'],4:['H','N','Q','W'],5:['R','K'],6:['D','E'],7:['C']}


# In[4]:


def _Str2Num(proteinsequence):
    """
    translate the amino acid letter into the corresponding class based on the
    
    given form.
    
    """
    repmat={}
    for i in _repmat:
        for j in _repmat[i]:
            repmat[j]=i
            
    res=proteinsequence
    for i in repmat:
        res=res.replace(i,str(repmat[i]))
    return res


# In[5]:


def CalculateConjointTriad(proteinsequence):
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    _repmat={1:["A",'G','V'],2:['I','L','F','P'],3:['Y','M','T','S'],4:['H','N','Q','W'],5:['R','K'],6:['D','E'],7:['C']}
    res={}
    proteinnum=_Str2Num(proteinsequence)
    for i in range(8):
        for j in range(8):
            for k in range(8):
                temp=str(i)+str(j)+str(k)
                res[temp]=proteinnum.count(temp)
    return res


# In[ ]:




