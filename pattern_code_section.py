import os
import re
import argparse
import shutil
from API import *
from project import *

def extract_code_from_response(generated_text):
    """从API响应中提取代码块"""
    pattern = r"```c\s*(.*?)\s*```"
    matches = re.findall(pattern, generated_text, re.DOTALL)
    if matches:
        return matches[0]
    return generated_text  # 如果没有代码块标记，则返回整个文本

def optimize_function(function_name, source_folder, output_folder, system_prompt_base, additional_prompt=""):
    """优化单个函数，添加OpenMP并行原语"""
    # 读取函数源码
    try:
        with open(f'{source_folder}/{function_name}.c', 'r') as f:
            function_content = f.read()
    except FileNotFoundError:
        print(f"警告: 找不到函数 {function_name} 的源文件，跳过优化")
        return False
    
    # 构建系统提示词
    system_prompt = system_prompt_base
    if additional_prompt:
        system_prompt += "\n" + additional_prompt
    
    # 调用API生成优化后的代码
    generated_text = sys_question(system_prompt, function_content)
    optimized_code = extract_code_from_response(generated_text)
    
    # 确保输出目录存在
    os.makedirs(output_folder, exist_ok=True)
    
    # 保存优化后的代码
    with open(f'{output_folder}/{function_name}.c', 'w') as f:
        f.write(optimized_code)
    
    print(f"函数 {function_name} 已优化并保存")
    return True

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='添加OpenMP并行原语到NPB基准测试')
    parser.add_argument(
        '--folder', 
        nargs='+', 
        default=['BT'],
        help='List of folder names to process (e.g., BT CG SP)'
    )
    parser.add_argument(
        '--class', 
        default='W',
        help='Problem class size (S, W, A, B, C)'
    )
    
    args = parser.parse_args()
    folders = args.folder
    class_size = getattr(args, 'class')  # 因为'class'是Python关键字
    
    base_dir = os.environ.get("BASE_DIR")
    npb_base_folder = f'{base_dir}/NPB3.0-omp-C'
    
    # 读取基础系统提示词
    with open(f"{base_dir}/prompt/0_code.txt", 'r') as f:
        system_prompt_base = f.read()
    
    # 遍历每个文件夹进行优化
    for folder_name in folders:
        bench_lower = folder_name.lower()
        source_folder = f'{npb_base_folder}/{folder_name}/function_#_omp'
        baseline_folder = f'{npb_base_folder}/{folder_name}/function_pattern_baseline'
        output_folder = f'{baseline_folder}/C_code'
        
        # 确保输出目录存在
        os.makedirs(output_folder, exist_ok=True)
        
        print(f"开始处理文件夹: {folder_name}")
        
        # 分析C文件中的函数调用关系
        c_file_path = f"{npb_base_folder}/{folder_name}/{bench_lower}_#_omp.c"
        all_functions, call_dependencies = analyze_c_function_calls(c_file_path)
        
        # 初始化NPB并获取串行执行时间
        init_NPB(folder_name)
        success, init_time = run_NPB(bench_lower, class_size)
        
        if not success:
            print(f"警告: {folder_name} 初始化运行失败，跳过该文件夹")
            continue
        
        print(f"{folder_name} 初始串行执行时间: {init_time}秒")
        best_time = init_time
        
        # 找出没有依赖的函数
        independent_functions = single_function(all_functions, call_dependencies)
        print(f"找到的没有依赖的函数: {independent_functions}")
        
        # 优化没有依赖的函数
        for function_name in independent_functions:
            print(f"优化没有依赖的函数: {function_name}")
            
            # 尝试优化，最多重试2次
            for attempt in range(3):  # 最多3次尝试（1次原始+2次重试）
                optimize_function(function_name, source_folder, output_folder, system_prompt_base)
                
                # 替换函数实现并运行
                replace_NPB(folder_name, function_name, f'function_pattern_baseline/C_code')
                success, current_time = run_NPB(bench_lower, class_size)
                
                if success and current_time < best_time:
                    print(f"函数 {function_name} 优化成功! 时间从 {best_time}秒 减少到 {current_time}秒")
                    best_time = current_time
                    break  # 成功优化，跳出重试循环
                else:
                    if attempt < 2:  # 如果还有重试机会
                        print(f"函数 {function_name} 优化尝试 {attempt+1} 失败，重试...")
                    else:
                        print(f"函数 {function_name} 优化失败，恢复原始实现")
                        # 恢复原始实现
                        replace_NPB(folder_name, function_name, 'function_#_omp')
        
        # 获取BFS顺序并逆序遍历
        bfs_order = get_bfs_order(all_functions, call_dependencies)
        bfs_order.reverse()
        
        # 优化有依赖的函数
        for function_name in bfs_order:
            # 检查该函数是否有依赖函数
            if function_name in call_dependencies and call_dependencies[function_name]:
                dependent_functions = call_dependencies[function_name]
                print(f"处理有依赖的函数: {function_name}，依赖函数: {dependent_functions}")
                
                # 保存当前状态（用于失败回滚）
                success, current_best_time = run_NPB(bench_lower, class_size)
                if not success:
                    current_best_time = best_time
                
                # 追踪成功添加并行的依赖函数
                successfully_parallelized = []
                
                # 首先优化依赖函数
                for dep_func in dependent_functions:
                    print(f"优化依赖函数: {dep_func}")
                    additional_prompt = "The outer layer has already used \"parallel\""
                    
                    if optimize_function(dep_func, source_folder, output_folder, system_prompt_base, additional_prompt):
                        replace_NPB(folder_name, dep_func, f'function_pattern_baseline/C_code')
                        success, time_after_dep = run_NPB(bench_lower, class_size)
                        
                        if success and time_after_dep < 9999:  # 运行成功就保留
                            print(f"依赖函数 {dep_func} 优化成功!")
                            successfully_parallelized.append(dep_func)
                        else:
                            print(f"依赖函数 {dep_func} 优化失败，恢复原始实现")
                            replace_NPB(folder_name, dep_func, 'function_#_omp')
                
                # 然后优化主函数
                print(f"优化主函数: {function_name}")
                additional_prompt = "The following functions need to add \"parallel\": " + ", ".join(successfully_parallelized) if successfully_parallelized else ""
                
                if optimize_function(function_name, source_folder, output_folder, system_prompt_base, additional_prompt):
                    replace_NPB(folder_name, function_name, f'function_pattern_baseline/C_code')
                    success, time_after_main = run_NPB(bench_lower, class_size)
                    
                    if success and time_after_main < current_best_time:
                        print(f"主函数 {function_name} 及其依赖优化成功! 时间从 {current_best_time}秒 减少到 {time_after_main}秒")
                        best_time = time_after_main
                    else:
                        print(f"主函数 {function_name} 优化失败，恢复原始实现")
                        # 恢复主函数及其所有依赖函数
                        replace_NPB(folder_name, function_name, 'function_#_omp')
                        for dep_func in successfully_parallelized:
                            replace_NPB(folder_name, dep_func, 'function_#_omp')
        
        print(f"{folder_name} 优化完成，最终执行时间: {best_time}秒")

if __name__ == "__main__":
    main()