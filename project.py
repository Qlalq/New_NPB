import os
import shutil
import subprocess
import re
import json
from collections import deque
base_dir = os.environ.get("BASE_DIR")
NPB_dir = f"{base_dir}/NPB3.0-omp-C"


def init_NPB(bench, mode="no_omp"):
    """
    NPB bench初始化
    
    Args:
        bench (str): NPB 基准测试名称，例如 "BT"
        mode (str): 初始化模式
            - "no_omp": 使用 {bench_lower}_#_omp.c (默认)
            - "ori": 使用 {bench_lower}_ori.c  
            - "autoPar": 使用 rose_{bench_lower}_#_omp.c
    """
    bench_lower = bench.lower()
    
    # 根据模式确定源文件路径
    if mode == "no_omp":
        src_file = f"{NPB_dir}/{bench}/{bench_lower}_#_omp.c"
    elif mode == "ori":
        src_file = f"{NPB_dir}/{bench}/{bench_lower}_ori.c"
    elif mode == "autoPar":
        src_file = f"{NPB_dir}/{bench}/rose_{bench_lower}_#_omp.c"
    elif mode == "aaai":
        src_file = f"{NPB_dir}/{bench}/{bench_lower}_aaai.c"
    else:
        raise ValueError(f"不支持的模式: {mode}.")
    
    dst_file = f"{NPB_dir}/{bench}/{bench_lower}.c"
    shutil.copy(src_file, dst_file)
    
    print(f"{bench} has been initialized with mode '{mode}'")

def replace_NPB(bench, function, folder):
    """
    替换 NPB 基准测试中的函数实现
    
    Args:
        bench (str): NPB 基准测试名称，例如 "BT"
        function (str): 需要替换的函数名称，例如 "add"
        folder (str): 函数所在的子目录，例如 "function_refinement"
    """
    import os
    
    base_dir = os.environ.get("BASE_DIR")
    NPB_dir = f"{base_dir}/NPB3.0-omp-C"
    bench_lower = bench.lower()
    
    src_file = f"{NPB_dir}/{bench}/{folder}/{function}.c"
    match_file = f"{NPB_dir}/{bench}/function_#_omp/{function}.c"
    dst_file = f"{NPB_dir}/{bench}/{bench_lower}.c"
    
    # 检查源文件是否存在
    if not os.path.exists(src_file):
        print(f"源文件不存在: {src_file}，跳过替换")
        return
    
    # 检查匹配文件是否存在
    if not os.path.exists(match_file):
        print(f"匹配文件不存在: {match_file}，跳过替换")
        return
    
    # 检查目标文件是否存在
    if not os.path.exists(dst_file):
        print(f"目标文件不存在: {dst_file}，跳过替换")
        return
    
    try:
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
    except Exception as e:
        print(f"替换过程中出错: {e}，跳过替换")


def run_NPB(bench, CLASS, timeout=300, times=1):
    """
    执行 NPB , 返回结果正确性 & 执行时间
    
    Args:
        bench (str): NPB 基准测试名称，例如 "BT"
        CLASS (str): 数据大小，例如 "S", 从小到大可选"S" "W" "A" "B" "C" 
        timeout (int): 执行超时时间(秒)，默认300秒
        times (int): 重复执行次数，默认1次
    
    Return:
        正常执行?/执行时间(s), 执行失败则执行时间是9999
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
    
    # 存储每次运行的时间
    execution_times = []
    
    # 重复执行 times 次
    for i in range(times):
        # 执行基准测试并获取输出
        try:
            result = subprocess.run(run_command, shell=True, check=True, cwd=NPB_dir, 
                                   capture_output=True, text=True, timeout=timeout)
            output = result.stdout
        except subprocess.TimeoutExpired:
            return False, 9999  # 执行超时，返回错误
        except subprocess.CalledProcessError:
            return False, 9999  # 运行失败，返回错误

        # 检查输出以确认结果
        verification_successful = re.search(r"Verification\s*=\s*SUCCESSFUL", output)
        time_match = re.search(r"Time in seconds\s*=\s*([\d.]+)", output)
        
        if verification_successful and time_match:
            time_taken = float(time_match.group(1))
            execution_times.append(time_taken)
        else:
            return False, 9999  # 验证失败，返回错误
    
    # 计算平均执行时间
    average_time = sum(execution_times) / len(execution_times)
    
    # 汇总打印结果
    print(f"{bench}.{CLASS} executed {times} times:")
    print(f"Times: {execution_times}")
    print(f"Average: {average_time:.3f}s")
    
    return True, average_time



def code_line(address):
    """
    返回文件中每一行的行号和内容。
    
    Args:
    address (str): 文件的路径。
    
    Return:
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



def analyze_c_function_calls(c_file_path):
    """
    分析C文件中自定义函数的包含关系。

    Args:
        c_file_path (str): C文件的路径。

    Return:
        tuple: 包含两个元素的元组:
               1. 一个集合，包含所有自定义函数的名称
               2. 一个字典，键是调用者函数名，值是被调用函数名组成的列表。
                  例如: {'error_norm': ['exact_solution']}
    """
    try:
        with open(c_file_path, 'r', encoding='utf-8') as f:
            code = f.read()
    except FileNotFoundError:
        print(f"错误：文件 {c_file_path} 未找到。")
        return set(), {}
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return set(), {}

    # 1. 找出所有自定义函数的名称
    #    正则表达式：匹配函数定义，如 "void func_name(...)" 或 "int func_name(...)"
    #    它会尝试匹配：
    #    - 可选的 static
    #    - 返回类型 (单词，可以带空格和*)
    #    - 函数名 (我们要捕获这个)
    #    - 括号和里面的参数 (我们不关心具体参数)
    #    - 必须以 { 结尾（表示函数体的开始）
    #    我们用 re.MULTILINE 让 ^ 匹配每行的开始
    function_definition_pattern = re.compile(
        r"^\s*(?:static\s+)?[\w\s\*]+\s+([a-zA-Z_]\w*)\s*\([^)]*\)\s*\{",
        re.MULTILINE
    )
    
    defined_function_names = set()
    # 使用 finditer 来获取匹配对象，从中提取函数名
    for match in function_definition_pattern.finditer(code):
        func_name = match.group(1)
        # 避免将 C 语言的关键字（如 if, for, while）误认为函数名
        # 实际中这个列表可能需要更完善
        if func_name not in ["if", "for", "while", "switch", "return", "sizeof"]: 
            defined_function_names.add(func_name)

    if not defined_function_names:
        print("在文件中没有找到函数定义。")
        return set(), {}
    
    print(f"找到的待优化函数有(含main): {defined_function_names}")

    call_graph = {} # 用来存放结果，比如 {'大厨A': ['小助手B', '小助手C']}

    # 2. 遍历每个找到的自定义函数，分析其函数体
    for caller_match in function_definition_pattern.finditer(code):
        caller_name = caller_match.group(1)

        # 我们只关心我们自己定义的函数之间的调用
        if caller_name not in defined_function_names and caller_name != "main": # main 也需要分析
            continue

        # 获取函数体的开始位置（即 { 之后）
        body_start_index = caller_match.end()
        
        # 寻找函数体的结束位置 (对应的 '}')
        # 这是一个简化的方法，通过匹配花括号对来确定函数体范围
        open_braces = 1
        current_pos = body_start_index
        body_end_index = -1
        while current_pos < len(code) and open_braces > 0:
            if code[current_pos] == '{':
                open_braces += 1
            elif code[current_pos] == '}':
                open_braces -= 1
            current_pos += 1
        
        if open_braces == 0:
            body_end_index = current_pos -1 # `current_pos` 现在在 `}` 的下一个字符
        else:
            # 如果括号不匹配，说明C代码可能有问题，或者我们的匹配逻辑不够完善
            print(f"警告: 函数 {caller_name} 的括号可能不匹配，分析可能不准确。")
            continue # 跳过这个函数

        function_body = code[body_start_index:body_end_index]

        # 3. 在函数体内查找对其他已定义函数的调用
        #    正则表达式：匹配函数调用，如 "func_name("
        #    \b 是单词边界，确保我们匹配到的是完整的函数名，而不是某个变量名的一部分
        function_call_pattern = re.compile(r"\b([a-zA-Z_]\w*)\s*\(")
        
        called_functions_in_body = set()
        for call_match in function_call_pattern.finditer(function_body):
            callee_name = call_match.group(1)
            # 如果被调用的函数是我们定义过的，并且不是调用自己（递归调用可以另外考虑）
            if callee_name in defined_function_names and callee_name != caller_name:
                called_functions_in_body.add(callee_name)
        
        if called_functions_in_body:
            if caller_name not in call_graph:
                call_graph[caller_name] = []
            call_graph[caller_name].extend(list(called_functions_in_body))
            # 去重，以防同一个函数在不同地方被多次添加到列表
            call_graph[caller_name] = sorted(list(set(call_graph[caller_name])))

    # 返回两个值：defined_function_names 和 call_graph
    return defined_function_names, call_graph



def get_bfs_order(all_functions, call_dependencies):
    """
    基于 def analyze_c_function_calls(c_file_path) 的优化顺序
    Args:
        all_functions (set or list): 所有自定义函数的名称集合或列表。
                                     例如: {'main', 'helper', 'util'}
        call_dependencies (dict): 函数调用关系的字典，
                                 键是调用者函数名，值是被调用函数名组成的列表。
                                 例如: {'main': ['helper'], 'helper': ['util']}
    Returns:
        list: 一个表示BFS遍历顺序的函数名列表。
              这个顺序满足拓扑排序的性质（前面的函数不能被后面的函数调用）。
    """
    
    graph = {func: call_dependencies.get(func, []) for func in all_functions}
    
    # 找到入度为0的节点作为起始节点(通常是main)
    # 如果没有显式指定main，则选择第一个没有被其他函数调用的函数
    in_degree = {func: 0 for func in all_functions}
    for caller, callees in call_dependencies.items():
        for callee in callees:
            if callee in in_degree:
                in_degree[callee] += 1
    
    # 如果main在函数集合中，从main开始BFS
    start_node = 'main' if 'main' in all_functions else None
    
    # 如果没有main，则选择第一个入度为0的函数
    if start_node is None:
        for func, degree in in_degree.items():
            if degree == 0:
                start_node = func
                break
    
    # 如果仍然没有找到起始节点，使用任意函数
    if start_node is None and all_functions:
        start_node = next(iter(all_functions))
    
    # 执行BFS
    result = []
    visited = set()
    queue = deque([start_node])
    visited.add(start_node)
    
    while queue:
        current = queue.popleft()
        result.append(current)
        
        # 获取当前函数调用的所有函数
        callees = graph.get(current, [])
        for callee in callees:
            if callee in all_functions and callee not in visited:
                queue.append(callee)
                visited.add(callee)
    
    # 添加剩余未访问的函数
    for func in all_functions:
        if func not in visited:
            result.append(func)
    
    return result



def test_function_Importance(Analysis_type='ori', folders=['BT'], CLASS='S'):
    """
    测试函数重要性分析并输出结果到 JSONL 文件。

    Args:
        Analysis_type (str): 分析文件夹后缀，默认是 'ori'。
        folders (list): 文件夹名称列表，默认是 ['BT']。
        CLASS (str): 数据大小，默认是 'S'。
    """
    base_dir = os.environ.get("BASE_DIR")
    NPB_dir = f"{base_dir}/NPB3.0-omp-C"
    output_file_path = f"{NPB_dir}/All/function_Importance_{Analysis_type}.jsonl"
    
    # 确保输出文件夹存在
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    results = []

    for folder in folders:
        # 1. 初始化 NPB
        init_NPB(folder)

        # 2. 分析 C 文件中的函数调用关系
        c_file_path = f"{NPB_dir}/{folder}/{folder.lower()}_#_omp.c"
        defined_functions, call_graph = analyze_c_function_calls(c_file_path)

        # 3. 获取 BFS 遍历顺序
        bfs_order = get_bfs_order(defined_functions, call_graph)

        # 反转 BFS 遍历顺序
        bfs_order.reverse()  # 或者使用 bfs_order = bfs_order[::-1]

        previous_time = None  # 保存上次运行的时间

        # 4. 遍历函数顺序并替换和运行 NPB
        for function in bfs_order:
            # 替换函数实现
            replace_NPB(folder, function, f'function_{Analysis_type}')

            # 运行 NPB 并获取结果
            success, execution_time = run_NPB(folder, CLASS)

            # 生成结果字典
            saved_time = 0 if previous_time is None else previous_time - execution_time
            result_entry = {
                "bench": folder,
                "function": function,
                "class": CLASS,
                "time": execution_time,
                "saved_time": saved_time
            }
            results.append(result_entry)
            print('*' * 10)
            print(f"bench:{folder}\nfunction:{function}\nclass:{CLASS}\ntime:{execution_time}\nsaved_time:{saved_time}")
            print('*' * 10)
            
            # 更新上次运行的时间
            if success:
                previous_time = execution_time
            else:
                previous_time = None  # 运行失败，重置为 None

    # 将结果写入 JSONL 文件
    with open(output_file_path, 'w') as jsonl_file:
        for result in results:
            jsonl_file.write(json.dumps(result) + '\n')

    print(f"结果已写入 {output_file_path}")


def single_function(all_functions, call_relations):
    """
    找出没有依赖的函数（不被其他函数调用的函数）。
    自动剔除main相关的数据。
    
    Args:
        all_functions (set): 包含所有自定义函数的名称的集合
        call_relations (dict): 键是调用者函数名，值是被调用函数名组成的列表
        
    Returns:
        set: 没有依赖的函数集合（不被其他函数调用的函数）
    """
    # 剔除main相关的数据
    if 'main' in all_functions:
        all_functions.remove('main')
    if 'main' in call_relations:
        del call_relations['main']
    
    # 收集所有被调用的函数，需要考虑keys和values
    called_functions = set()
    
    # 添加所有被调用的函数(values)
    for called_list in call_relations.values():
        called_functions.update(called_list)
    
    # 添加所有调用其他函数的函数(keys)
    called_functions.update(call_relations.keys())
    
    # 计算没有依赖的函数（不被其他函数调用且不调用其他函数的函数）
    independent_functions = all_functions - called_functions
    
    return independent_functions


