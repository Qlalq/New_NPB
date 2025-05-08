import os
base_dir = os.environ.get("BASE_DIR")

# 定义源文件路径和目标文件夹路径
source_file = f'{base_dir}/NPB3.0-omp-C/BT/bt.c'
target_folder = f'{base_dir}/NPB3.0-omp-C/BT/function_#_omp'


# 确保目标文件夹存在
os.makedirs(target_folder, exist_ok=True)

# 定义要提取的函数列表
functions = [
    'static void add(void)',
    'static void adi(void)',
    'static void error_norm(double rms[5])',
    'static void rhs_norm(double rms[5])',
    'static void exact_rhs(void)',
    'static void exact_solution(double xi, double eta, double zeta, double dtemp[5])',
    'static void initialize(void)',
    'static void lhsinit(void)',
    'static void lhsx(void)',
    'static void lhsy(void)',
    'static void lhsz(void)',
    'static void compute_rhs(void)',
    'static void set_constants(void)',
    'static void verify(int no_time_steps, char *class, boolean *verified)',
    'static void x_solve(void)',
    'static void x_backsubstitute(void)',
    'static void x_solve_cell(void)',
    'static void matvec_sub(double ablock[5][5], double avec[5], double bvec[5])',
    'static void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5])',
    'static void binvcrhs(double lhs[5][5], double c[5][5], double r[5])',
    'static void binvrhs(double lhs[5][5], double r[5])',
    'static void y_solve(void)',
    'static void y_backsubstitute(void)',
    'static void y_solve_cell(void)',
    'static void z_solve(void)',
    'static void z_backsubstitute(void)',
    'static void z_solve_cell(void)'
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
    
    # 找到函数体的开始和结束位置
    start_index = content.find('{', last_index)
    if start_index == -1:
        print(f"Function body start not found for {func}.")
        continue
    
    # 使用栈来匹配大括号
    stack = []
    end_index = start_index
    for i in range(start_index, len(content)):
        if content[i] == '{':
            stack.append('{')
        elif content[i] == '}':
            stack.pop()
            if not stack:
                end_index = i
                break
    
    if not stack:
        # 提取函数体内容
        function_body = content[last_index:end_index + 1]
        
        # 保存到文件
        func_name = func.split('(')[0].split()[-1]
        target_file = os.path.join(target_folder, f"{func_name}.c")
        with open(target_file, 'w') as f:
            f.write(function_body)
        print(f"Function {func_name} saved to {target_file}")
    else:
        print(f"Function body end not found for {func}.")