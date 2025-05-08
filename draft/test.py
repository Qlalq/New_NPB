import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from API import ask_gpt_question, patch_sys_question
from project import code_line

base_dir = os.environ.get("BASE_DIR")
NPB_dir = f"{base_dir}/NPB3.0-omp-C"
file_dir = f'{base_dir}/NPB3.0-omp-C/BT/function_#_omp/compute_rhs.c'

print(code_line(file_dir))
