#%%
import numpy as np
import csv

data20070618 = np.load('20070618.npz')
ice_area20070618 = data20070618['icebergArea']
survey_area20070618 = data20070618['surveyArea']

data20070822 = np.load('20070822.npz')
ice_area20070822 = data20070822['icebergArea']
survey_area20070822 = data20070822['surveyArea']

print(sum(ice_area20070618)/survey_area20070618)
print(sum(ice_area20070822)/survey_area20070822)

print((ice_area20070618 > 1.6).sum())
print((ice_area20070822 > 1.6).sum())

np.histogram(ice_area20070822)

# %%
