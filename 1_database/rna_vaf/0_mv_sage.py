#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import os
import pandas as pd
from helpers import *


# #### Go!

# In[13]:


I_FILES = [i for i in os.popen('gsutil ls gs://hmf-rna-analysis/sage_append/somatic/').read().split("\n") if 'somatic' in i]


# In[14]:


for f in I_FILES:
    cmd = "gsutil cp " + f + " " + SAGE_DIR
    if not os.path.exists(SAGE_DIR + f.split("/somatic/")[1]):
        os.system(cmd)

