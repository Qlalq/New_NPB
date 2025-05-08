"""
测试文件夹内函数正确性
"""

from API import sys_question, patch_sys_question
from project import init_NPB, replace_NPB, run_NPB
import os
import json

base_dir = os.environ.get("BASE_DIR")
file_address = f"{base_dir}/NPB3.0-omp-C/BT/function_refinement/compute_rhs.c"


NPB_dir = f"{base_dir}/NPB3.0-omp-C"
result_file = f"{NPB_dir}/All/refine.jsonl"


benches = ["BT", "CG", "EP", "FT", "LU", "MG", "SP"]

for bench in benches:
    bench_dir = f"{NPB_dir}/{bench}/function_refinement"
    for function in os.listdir(bench_dir):
        init_NPB(bench)
        if function.endswith(".c"):
            function_name = os.path.splitext(function)[0]
            replace_NPB(bench, function_name, "function_refinement")
            work, time = run_NPB(bench, "S")
            print('*'*10)
            print(work, time)
            print('*'*10)
            result = {
                "bench": bench,
                "function": function_name,
                "CLASS": "S",
                "work": work,
                "time": time
            }
            with open(result_file, "a") as f:
                json.dump(result, f)
                f.write("\n")
            print(f"Processed {bench}/{function_name}")