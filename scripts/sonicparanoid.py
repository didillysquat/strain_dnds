# output argument will always be sonicparanoid/output
# input will always be sonicparanoid

import os

# first check to see if there is already a run called parkinson in the sonicparanoid/output folder
print("the current working directory is os.getcwd()")
print("the input directory is " + os.path.abspath("sonicparanoid"))
print("the output directory is " + os.path.abspath("sonicparanoid/output"))