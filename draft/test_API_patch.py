import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from API import *
from project import *

base_dir = os.environ.get("BASE_DIR")
before_address = f"{base_dir}/draft/add.c"
sys_address = f"{base_dir}/prompt/0_code.txt"
after_address = f"{base_dir}/draft/after.c"
with open(sys_address, 'r') as f:
    sys_prompt = f.read()
with open(before_address, 'r') as f:
    prompt = f.read()

result = sys_question(sys_prompt, prompt)
print(result)
with open(after_address, 'w') as f:
    f.write(result)
# prompt = code_line(before_address)
# # print(prompt)
# print(patch_sys_question(sys_prompt, prompt, before_address))