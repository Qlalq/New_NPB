import os
base_dir = os.environ.get("BASE_DIR")

# 定义源文件路径和目标文件夹路径
source_file = f'{base_dir}/NPB3.0-omp-C/CG/cg_ori.c'
target_folder = f'{base_dir}/NPB3.0-omp-C/CG/function_ori'

# 确保目标文件夹存在
os.makedirs(target_folder, exist_ok=True)

# 定义要提取的函数列表
functions = [
    'static void conj_grad',
    'static void makea',
    'static void sparse',
    'static void sprnvc',
    'static int icnvrt',
    'static void vecset'
]

# 读取源文件内容
with open(source_file, 'r') as file:
    content = file.read()

# 遍历每个函数并提取内容
for func in functions:
    # 找到最后一个函数定义的位置
    last_index = content.rfind(func)
    if last_index == -1:
        print(f"Function {func} not found in the source file.")
        continue
    
    # 找到参数部分的开始和结束位置
    param_start = content.find('(', last_index)
    if param_start == -1:
        print(f"Parameter start not found for {func}.")
        continue
    
    # 匹配参数部分的括号
    stack = []
    param_end = param_start
    for i in range(param_start, len(content)):
        if content[i] == '(':
            stack.append('(')
        elif content[i] == ')':
            stack.pop()
            if not stack:
                param_end = i
                break
    
    if not stack:
        # 提取参数部分
        function_header = content[last_index:param_end + 1]
        
        # 找到函数体的开始和结束位置
        body_start = content.find('{', param_end)
        if body_start == -1:
            print(f"Function body start not found for {func}.")
            continue
        
        # 匹配函数体的大括号
        stack = []
        body_end = body_start
        for i in range(body_start, len(content)):
            if content[i] == '{':
                stack.append('{')
            elif content[i] == '}':
                stack.pop()
                if not stack:
                    body_end = i
                    break
        
        if not stack:
            # 提取函数体内容
            function_body = content[last_index:body_end + 1]
            
            # 保存到文件
            func_name = func.split()[-1]  # 提取函数名
            target_file = os.path.join(target_folder, f"{func_name}.c")
            with open(target_file, 'w') as f:
                f.write(function_body)
            print(f"Function {func_name} saved to {target_file}")
        else:
            print(f"Function body end not found for {func}.")
    else:
        print(f"Parameter end not found for {func}.")