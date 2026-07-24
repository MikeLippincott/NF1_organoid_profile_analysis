#!/usr/bin/env python
# coding: utf-8

# This notebook unzips the data compressed file and sets up the data directory for analysis.
# This way all users have the same data directory structure and files for analysis.

# In[ ]:


import pathlib
import shutil

from notebook_init_utils import init_notebook

root_dir, in_notebook = init_notebook()


# In[7]:


OVERWRITE = False
shippable_zip_path = pathlib.Path(f"{root_dir}/data/shippable_dir.zip").resolve(
    strict=True
)


# In[9]:


unpack_dir = pathlib.Path(f"{root_dir}/data/").resolve(strict=False)
# check if unpack_dir is empty or not
if unpack_dir.exists() and any(unpack_dir.iterdir()) and not OVERWRITE:
    print("Data directory is not empty and OVERWRITE is False. Skipping unpacking.")
else:
    # unzip the shippable_dir.zip file to the data directory
    shutil.unpack_archive(shippable_zip_path, unpack_dir, "zip")
