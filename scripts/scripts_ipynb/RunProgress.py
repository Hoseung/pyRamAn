
# coding: utf-8

# #### Zoom in simulation progress (~40 clusters)

# ##### 1. retrieve most recent log files. 
# Can it be done in one go? 
# 
# 

# In[ ]:




# ##### 2. parse the log to find (last aexp, # of stellar particles, # of grids, memory usage)

# In[ ]:




# ##### 3. plot progress bars

# In[1]:

# data
clusters = ("01605", "29195")
aexp_not = [0.8, 0.4]


#%%
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()  # For second x ticks

# In[2]:
y_pos = np.arange(len(clusters))
ax1.barh(y_pos, aexp_not, align='center', alpha=0.4)
plt.yticks(y_pos, clusters)
#ax1.set_yticks(y_pos, clusters)
ax1.set_xlabel("exponential factor")
#%%
zreds = np.asarray([40, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.0])
aexps = 1 / (1+zreds)

ax2.set_xticks(aexps)
ax2.set_xticklabels(zreds)
ax2.set_xlabel("Redshift")

fig.show()
