import os
import shutil
import subprocess
import re

base_dir = os.environ.get("BASE_DIR")
NPB_dir = f"{base_dir}/NPB3.0-omp-C"


def init_NPB(bench):
    """
    NPB bench初始化
    
    Args:
        bench (str): NPB 基准测试名称，例如 "BT"
    """
    bench_lower = bench.lower()

    src_file = f"{NPB_dir}/{bench}/{bench_lower}_#_omp.c"
    dst_file = f"{NPB_dir}/{bench}/{bench_lower}.c"
    shutil.copy(src_file, dst_file)
    
    print(f"{bench} has been initialized")



def replace_NPB(bench, function, folder):
    """
    替换 NPB 基准测试中的函数实现
    
    Args:
        bench (str): NPB 基准测试名称，例如 "BT"
        function (str): 需要替换的函数名称，例如 "add"
        folder (str): 函数所在的子目录，例如 "function_refinement"
    """
    base_dir = os.environ.get("BASE_DIR")
    NPB_dir = f"{base_dir}/NPB3.0-omp-C"
    bench_lower = bench.lower()
    
    src_file = f"{NPB_dir}/{bench}/{folder}/{function}.c"
    match_file = f"{NPB_dir}/{bench}/function_#_omp/{function}.c"
    dst_file = f"{NPB_dir}/{bench}/{bench_lower}.c"
    
    with open(match_file, "r") as f:
        match_content = f.read()
    
    with open(src_file, "r") as f:
        new_content = f.read().strip()

    with open(dst_file, "r") as f:
        content = f.read()
    
    content = content.replace(match_content, new_content)
    
    with open(dst_file, "w") as f:
        f.write(content)
    
    print(f"已将 {function} 函数从 {folder} 替换到 {bench_lower}.c 中")




def run_NPB(bench, CLASS):
    """
    执行 NPB , 返回结果正确性 & 执行时间
    
    Args:
        bench (str): NPB 基准测试名称，例如 "BT"
        CLASS (str): 数据大小，例如 "S", 从小到大可选"S" "W" "A" "B" "C" 
    """
    bench = bench.lower()
    # 构建 make 命令
    make_command = f"make {bench} CLASS={CLASS}"
    
    # 切换到 NPB 目录并执行 make 命令
    try:
        subprocess.run(make_command, shell=True, check=True, cwd=NPB_dir)
    except subprocess.CalledProcessError:
        return False, 9999  # make 失败，返回错误
    
    # 构建运行命令
    run_command = f"./bin/{bench}.{CLASS}"
    
    # 执行基准测试并获取输出
    try:
        result = subprocess.run(run_command, shell=True, check=True, cwd=NPB_dir, capture_output=True, text=True)
        output = result.stdout
    except subprocess.CalledProcessError:
        return False, 9999  # 运行失败，返回错误

    # 检查输出以确认结果
    verification_successful = re.search(r"Verification\s*=\s*SUCCESSFUL", output)
    time_match = re.search(r"Time in seconds\s*=\s*([\d.]+)", output)
    
    if verification_successful and time_match:
        time_taken = float(time_match.group(1))
        return True, time_taken  # 返回结果正确和执行时间
    else:
        return False, 9999  # 验证失败，返回错误


def code_line(address):
    """
    返回文件中每一行的行号和内容。
    
    参数:
    address (str): 文件的路径。
    
    返回:
    str: 包含每一行行号和内容的多行字符串。
    """
    try:
        with open(address, 'r') as file:
            lines = file.readlines()
            output = []
            for i, line in enumerate(lines, start=1):
                if line.strip():  # 检查是否为空行
                    output.append(f"{i}. {line.rstrip()}")
    except FileNotFoundError:
        return f"无法找到文件: {address}"
    except IOError:
        return f"读取文件 {address} 时出错"
    
    return '\n'.join(output)